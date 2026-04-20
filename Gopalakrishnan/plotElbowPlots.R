suppressPackageStartupMessages({
  library(ape)
  library(CATMicrobiome)
  library(coin)
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(mclust)
})

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  idx <- grep("^--file=", cmd_args)
  if (length(idx)) return(normalizePath(sub("^--file=", "", cmd_args[idx[1]])))
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(normalizePath(ofile))
  NA_character_
}

normalize_provider <- function(x) {
  match.arg(tolower(x), c("openai", "claude", "gemini", "combined"))
}

provider_dir_name <- function(provider) {
  switch(
    provider,
    openai = "openAIGenerated",
    claude = "claudeGenerated",
    gemini = "geminiGenerated",
    combined = "combined"
  )
}

parse_env_list <- function(name, default_values) {
  raw_value <- Sys.getenv(name, unset = "")
  if (!nzchar(raw_value)) {
    return(default_values)
  }

  values <- trimws(strsplit(tolower(raw_value), ",", fixed = TRUE)[[1]])
  values <- values[nzchar(values)]
  invalid_values <- setdiff(values, default_values)
  if (length(invalid_values)) {
    stop("Invalid values for ", name, ": ", paste(invalid_values, collapse = ", "))
  }

  values
}

parse_env_seed <- function(name, default_value = 123L) {
  raw_value <- Sys.getenv(name, unset = as.character(default_value))
  seed_value <- suppressWarnings(as.integer(raw_value))
  if (is.na(seed_value)) {
    stop("Invalid integer seed for ", name, ": ", raw_value)
  }
  seed_value
}

rank_display_name <- function(rank) {
  rank_norm <- tolower(rank)
  switch(rank_norm, otu = "OTU", tools::toTitleCase(rank_norm))
}

selection_target_folder <- function(target) {
  switch(
    tolower(target),
    family = "selectionTargetFamily",
    genus = "selectionTargetGenus",
    species = "selectionTargetSpecies",
    otu = "selectionTargetOTU",
    stop("Unsupported target: ", target)
  )
}

clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    str_replace("^.*?__", "") %>%
    str_replace_all("\\[|\\]", "") %>%
    str_squish()
}

script_path <- get_script_path()
script_dir <- if (!is.na(script_path)) dirname(script_path) else normalizePath(getwd())
provider <- normalize_provider(Sys.getenv("LLM_PROVIDER", unset = "combined"))
provider_dir <- provider_dir_name(provider)
run_seed <- parse_env_seed("LLM_SEED", 123L)
options(chai.seed = run_seed)
set.seed(run_seed)

output_root <- file.path(script_dir, "plots", "elbowPlot", provider)
cache_root <- Sys.getenv(
  "LLM_ELBOW_CACHE_ROOT",
  unset = file.path(script_dir, "cache", "elbowPlot", provider)
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)

source(file.path(script_dir, "code", "chai.R"))

build_ps <- function() {
  otu_path <- system.file("extdata", "d1OTUtable.csv", package = "CATMicrobiome")
  tax_path_local <- file.path(script_dir, "LLMCode", "d1Taxonomy.csv")
  tax_path <- if (file.exists(tax_path_local)) {
    tax_path_local
  } else {
    system.file("extdata", "d1Taxonomy.csv", package = "CATMicrobiome")
  }
  meta_path <- system.file("extdata", "d1Meta.csv", package = "CATMicrobiome")
  tree_path <- system.file("extdata", "d1Tree.tree", package = "CATMicrobiome")

  otutable <- read.csv(otu_path, header = TRUE, row.names = 1)
  taxonomy <- read.csv(tax_path, header = TRUE, row.names = 1)
  meta_data <- read.csv(meta_path, header = TRUE, row.names = 1)
  tree <- read.tree(tree_path)

  otus_keep <- intersect(tree$tip.label, rownames(taxonomy))
  tree <- keep.tip(tree, otus_keep)
  taxonomy <- taxonomy[match(tree$tip.label, rownames(taxonomy)), , drop = FALSE]
  otutable <- otutable[otus_keep, , drop = FALSE]

  phyloseq(
    otu_table(as.matrix(otutable), taxa_are_rows = TRUE),
    tax_table(as.matrix(taxonomy)),
    sample_data(meta_data),
    phy_tree(tree)
  )
}

build_rank_ps <- function(ps, rank) {
  rank_norm <- tolower(rank)
  if (rank_norm == "otu") return(ps)
  if (rank_norm == "species") {
    return(tax_glom(ps, taxrank = "species", NArm = FALSE))
  }
  tax_glom(ps, taxrank = rank_norm)
}

compute_pz_wilcox <- function(t_abun, outcome) {
  stopifnot(nrow(t_abun) == length(outcome))
  p_z_values <- sapply(colnames(t_abun), function(feature) {
    feature_values <- t_abun[, feature, drop = FALSE]
    wilcox <- wilcox_test(feature_values[, 1] ~ outcome)
    c(
      p_value = as.numeric(pvalue(wilcox)),
      z_value = as.numeric(statistic(wilcox))
    )
  }, simplify = "matrix")

  out <- as.data.frame(t(p_z_values))
  out$feature <- rownames(out)
  rownames(out) <- NULL
  out
}

make_main_df <- function(t_abun, pz_df, ps_obj, rank) {
  if (tolower(rank) == "otu") {
    pz_map <- if ("feature" %in% names(pz_df)) {
      pz_df %>%
        transmute(
          otu_rep = feature,
          p_main = as.numeric(p_value),
          z_main = as.numeric(z_value)
        )
    } else {
      data.frame(
        otu_rep = colnames(t_abun),
        p_main = as.numeric(pz_df$p_value),
        z_main = as.numeric(pz_df$z_value),
        stringsAsFactors = FALSE
      )
    }

    return(
      pz_map %>%
        mutate(
          rank_value = otu_rep,
          taxon_key = otu_rep
        )
    )
  }

  tax_df <- as.data.frame(tax_table(ps_obj)) %>%
    rownames_to_column("otu_rep")

  rank_col <- names(tax_df)[tolower(names(tax_df)) == tolower(rank)]
  if (length(rank_col) != 1) {
    stop("Could not uniquely find rank column '", rank, "' in tax_table.")
  }

  tax_map <- tax_df %>%
    transmute(otu_rep, rank_value = as.character(.data[[rank_col]]))

  pz_map <- if ("feature" %in% names(pz_df)) {
    pz_df %>%
      transmute(
        otu_rep = feature,
        p_main = as.numeric(p_value),
        z_main = as.numeric(z_value)
      )
  } else {
    data.frame(
      otu_rep = colnames(t_abun),
      p_main = as.numeric(pz_df$p_value),
      z_main = as.numeric(pz_df$z_value),
      stringsAsFactors = FALSE
    )
  }

  pz_map %>%
    left_join(tax_map, by = "otu_rep") %>%
    mutate(taxon_key = clean_taxon(rank_value))
}

build_target_main_df <- function(ps, rank) {
  rank_ps <- build_rank_ps(ps, rank)
  cs <- colSums(otu_table(rank_ps))
  cs[cs == 0] <- 1
  t_abun <- t(sweep(otu_table(rank_ps), 2, cs, "/") * 100)
  outcome <- factor(sample_data(rank_ps)$BinOutcomes, levels = c("NR", "R"))
  pz <- compute_pz_wilcox(t_abun, outcome)
  make_main_df(t_abun = t_abun, pz_df = pz, ps_obj = rank_ps, rank = rank) %>%
    select(taxon_key, p_main, z_main) %>%
    distinct()
}

prepare_aux_df <- function(aux_path, target_rank) {
  aux_df <- read_csv(aux_path, show_col_types = FALSE)

  key_col <- if (tolower(target_rank) == "otu") {
    names(aux_df)[tolower(names(aux_df)) == "otuname"]
  } else {
    names(aux_df)[tolower(names(aux_df)) == tolower(target_rank)]
  }

  if (length(key_col) != 1) {
    stop("Could not uniquely find key column for target rank '", target_rank, "' in ", aux_path)
  }

  aux_df %>%
    transmute(
      taxon_key = if (tolower(target_rank) == "otu") {
        as.character(.data[[key_col]])
      } else {
        clean_taxon(.data[[key_col]])
      },
      xA = as.numeric(A),
      xB = as.numeric(B),
      xC = as.numeric(C)
    ) %>%
    filter(!is.na(taxon_key), taxon_key != "") %>%
    distinct()
}

load_or_fit <- function(cache_file, fit_fun) {
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  obj <- fit_fun()
  saveRDS(obj, cache_file)
  obj
}

extract_bic_profile <- function(fit) {
  bic_mat <- as.matrix(fit$BIC)
  if (is.null(dim(bic_mat))) {
    stop("Mclust fit does not contain a BIC matrix.")
  }

  model_names <- colnames(bic_mat)
  profile_rows <- lapply(seq_len(nrow(bic_mat)), function(i) {
    row_vals <- bic_mat[i, ]
    if (all(is.na(row_vals))) {
      return(data.frame(
        G = as.integer(rownames(bic_mat)[i]),
        best_bic = NA_real_,
        best_model = NA_character_
      ))
    }
    best_idx <- which.max(row_vals)
    data.frame(
      G = as.integer(rownames(bic_mat)[i]),
      best_bic = as.numeric(row_vals[best_idx]),
      best_model = model_names[best_idx]
    )
  })

  bind_rows(profile_rows)
}

plot_elbow <- function(profile_df, selected_g, selected_model, title_text) {
  selected_row <- profile_df %>%
    filter(G == selected_g) %>%
    slice_head(n = 1)

  ggplot(profile_df, aes(x = G, y = best_bic)) +
    geom_line(color = "#2b6cb0", linewidth = 0.8, na.rm = TRUE) +
    geom_point(color = "#2b6cb0", size = 2.2, na.rm = TRUE) +
    geom_point(
      data = selected_row,
      color = "#c53030",
      size = 3.2,
      na.rm = TRUE
    ) +
    geom_text(
      data = selected_row,
      aes(label = paste0("Selected G=", selected_g, " (", selected_model, ")")),
      nudge_y = 0.03 * diff(range(profile_df$best_bic, na.rm = TRUE)),
      color = "#c53030",
      size = 3.5,
      fontface = "bold",
      na.rm = TRUE
    ) +
    scale_x_continuous(breaks = sort(unique(profile_df$G))) +
    labs(
      title = title_text,
      x = "Number of Gaussian Components (G)",
      y = "Best Mclust BIC"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

build_model_df <- function(main_df, aux_df) {
  main_df %>%
    left_join(aux_df, by = "taxon_key", relationship = "many-to-many") %>%
    mutate(across(c(xA, xB, xC), ~ tidyr::replace_na(.x, 0)))
}

provider_aux_file <- function(target_rank, aux_mode, aux_source = "phylum") {
  target_label <- if (tolower(target_rank) == "otu") "OTU" else tolower(target_rank)
  file.path(
    script_dir,
    "sideCov",
    provider_dir,
    aux_mode,
    selection_target_folder(target_rank),
    sprintf("d1Taxonomy_%s_ABCD_%s_and_below_%s.csv", target_label, aux_source, aux_mode)
  )
}

ps <- build_ps()
target_ranks <- parse_env_list("LLM_TARGET_RANKS", c("family", "genus", "species", "otu"))
aux_modes <- parse_env_list("LLM_AUX_MODES", c("weighted", "unweighted"))
aux_source <- Sys.getenv("LLM_ELBOW_AUX_SOURCE", unset = "phylum")
main_df_map <- setNames(lapply(target_ranks, function(rank) build_target_main_df(ps, rank)), target_ranks)

summary_rows <- vector("list", length(target_ranks) * length(aux_modes))
idx <- 1L

for (aux_mode in aux_modes) {
  out_dir <- file.path(output_root, aux_mode)
  mode_cache_dir <- file.path(cache_root, aux_mode)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mode_cache_dir, recursive = TRUE, showWarnings = FALSE)

  for (target_rank in target_ranks) {
    aux_file <- provider_aux_file(target_rank, aux_mode, aux_source)
    if (!file.exists(aux_file)) {
      warning("Skipping missing auxiliary file: ", aux_file)
      next
    }

    message("Generating elbow plot for mode=", aux_mode, ", target=", target_rank, ", aux=", aux_source)

    aux_df <- prepare_aux_df(aux_file, target_rank)
    df_model <- build_model_df(main_df_map[[target_rank]], aux_df)

    z <- df_model$z_main
    X <- as.data.frame(df_model[, c("xA", "xB", "xC"), drop = FALSE])
    keep <- complete.cases(z, X)
    z <- z[keep]
    X <- X[keep, , drop = FALSE]

    if (length(z) < 10) {
      warning("Skipping target rank ", target_rank, " / ", aux_mode, ": too few complete rows.")
      next
    }

    fit_df <- as.data.frame(X)
    names(fit_df) <- paste0("x", seq_len(ncol(fit_df)))
    fit_df$z <- z

    cache_file <- file.path(mode_cache_dir, sprintf("cache_mclust_target_%s_aux_%s.rds", target_rank, aux_source))
    mclust_fit <- load_or_fit(
      cache_file,
      function() {
        set.seed(run_seed)
        fit_mclust_with_fallback(fit_df, K_vec = 2:10)
      }
    )

    profile_df <- extract_bic_profile(mclust_fit)
    out_pdf <- file.path(
      out_dir,
      sprintf("elbow_target_%s_aux_%s_%s.pdf", target_rank, aux_source, aux_mode)
    )

    p <- plot_elbow(
      profile_df,
      selected_g = mclust_fit$G,
      selected_model = mclust_fit$modelName,
      title_text = paste0(
        "Gaussian Mixture Elbow Plot (",
        provider,
        ", ",
        aux_mode,
        ", target=",
        rank_display_name(target_rank),
        ", sideCov=",
        rank_display_name(aux_source),
        ")"
      )
    )
    ggsave(out_pdf, p, width = 7, height = 5)

    summary_rows[[idx]] <- tibble::tibble(
      provider = provider,
      aux_mode = aux_mode,
      aux_source = aux_source,
      target_rank = target_rank,
      selected_g = as.integer(mclust_fit$G),
      selected_model = as.character(mclust_fit$modelName),
      n_features = nrow(fit_df),
      output_pdf = normalizePath(out_pdf, mustWork = FALSE)
    )
    idx <- idx + 1L
  }
}

summary_df <- bind_rows(summary_rows)
summary_file <- file.path(output_root, "elbow_bic_summary.csv")
readr::write_csv(summary_df %>% arrange(aux_mode, target_rank), summary_file)
print(summary_df)
message("Saved elbow plot summary: ", summary_file)
