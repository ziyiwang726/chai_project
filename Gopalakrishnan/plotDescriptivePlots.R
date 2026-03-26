suppressPackageStartupMessages({
  library(ape)
  library(CATMicrobiome)
  library(coin)
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  idx <- grep("^--file=", cmd_args)
  if (length(idx)) return(normalizePath(sub("^--file=", "", cmd_args[idx[1]])))
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(normalizePath(ofile))
  NA_character_
}

script_path <- get_script_path()
script_dir <- if (!is.na(script_path)) dirname(script_path) else normalizePath(getwd())
input_root <- file.path(script_dir, "sideCov")
output_root <- file.path(script_dir, "plots", "descriptivePlots")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    str_replace("^.*?__", "") %>%
    str_replace_all("\\[|\\]", "") %>%
    str_squish()
}

rank_display_name <- function(rank) {
  rank_norm <- tolower(rank)
  switch(
    rank_norm,
    otu = "OTU",
    tools::toTitleCase(rank_norm)
  )
}

build_ps <- function() {
  otu_path <- system.file("extdata", "d1OTUtable.csv", package = "CATMicrobiome")
  tax_path <- system.file("extdata", "d1Taxonomy.csv", package = "CATMicrobiome")
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

parse_aux_metadata <- function(path) {
  rel <- sub(paste0("^", normalizePath(input_root), "/"), "", normalizePath(path))
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  if (length(parts) < 4) stop("Unexpected sideCov path: ", path)

  target_rank <- switch(
    parts[3],
    selectionTargetFamily = "family",
    selectionTargetGenus = "genus",
    selectionTargetSpecies = "species",
    selectionTargetOTU = "otu",
    stop("Unknown target folder: ", parts[3])
  )

  file_name <- basename(path)
  aux_source <- str_match(file_name, "_ABCD_([^_]+)_and_below_")[, 2]
  if (is.na(aux_source)) stop("Could not parse auxiliary source from ", file_name)

  list(
    provider_dir = parts[1],
    mode = parts[2],
    target_folder = parts[3],
    target_rank = target_rank,
    aux_source = aux_source,
    file_name = file_name,
    path = path
  )
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

marker_labels <- c(
  xA = "A: Melanoma response",
  xB = "B: Melanoma survival",
  xC = "C: Other cancer"
)

marker_colors <- c(
  xA = "#1b9e77",
  xB = "#d95f02",
  xC = "#7570b3"
)

plot_descriptive <- function(df_model, meta) {
  df_long <- df_model %>%
    select(taxon_key, z_main, xA, xB, xC) %>%
    pivot_longer(
      cols = c(xA, xB, xC),
      names_to = "panel",
      values_to = "aux_score"
    ) %>%
    filter(is.finite(z_main), is.finite(aux_score)) %>%
    group_by(panel) %>%
    filter(n() >= 3) %>%
    ungroup()

  if (!nrow(df_long)) return(NULL)

  x_lim <- max(abs(df_long$z_main), na.rm = TRUE)
  y_lim <- max(abs(df_long$aux_score), na.rm = TRUE)
  x_lim <- ceiling(x_lim * 10) / 10
  y_lim <- ceiling(y_lim * 10) / 10

  ggplot(df_long, aes(x = z_main, y = aux_score, color = panel)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_point(alpha = 0.68, size = 1.5) +
    geom_smooth(
      aes(fill = panel),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.75
    ) +
    scale_color_manual(values = marker_colors, labels = marker_labels) +
    scale_fill_manual(values = marker_colors, labels = marker_labels) +
    coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(-y_lim, y_lim), expand = FALSE) +
    labs(
      x = "Wilcoxon z-score",
      y = "Auxiliary z-score",
      color = "Marker"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 10, 8, 10)
    ) +
    guides(fill = "none")
}

aux_files <- list.files(
  input_root,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
aux_files <- aux_files[str_detect(aux_files, "selectionTarget")]

if (!length(aux_files)) {
  stop("No auxiliary files found under ", input_root)
}

ps <- build_ps()
target_ranks <- c("family", "genus", "species", "otu")
main_df_map <- setNames(lapply(target_ranks, function(rank) build_target_main_df(ps, rank)), target_ranks)

results <- vector("list", length(aux_files))

for (i in seq_along(aux_files)) {
  meta <- parse_aux_metadata(aux_files[i])
  message(
    "Processing mode=", meta$mode,
    ", target=", meta$target_rank,
    ", aux=", meta$aux_source
  )

  aux_df <- prepare_aux_df(meta$path, meta$target_rank)
  df_model <- main_df_map[[meta$target_rank]] %>%
    left_join(aux_df, by = "taxon_key", relationship = "many-to-many")

  plot_obj <- plot_descriptive(df_model, meta)
  if (is.null(plot_obj)) {
    warning("Skipping plot with insufficient complete rows for ", meta$file_name)
    next
  }

  out_dir <- file.path(
    output_root,
    meta$provider_dir,
    meta$mode,
    meta$target_folder
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_file <- file.path(
    out_dir,
    paste0(
      "wilcoxz_descriptive_",
      meta$target_rank,
      "_target_aux_",
      meta$aux_source,
      "_",
      meta$mode,
      ".pdf"
    )
  )

  ggsave(out_file, plot_obj, width = 7, height = 5)

  results[[i]] <- data.frame(
    mode = meta$mode,
    target_rank = meta$target_rank,
    aux_source = meta$aux_source,
    output = out_file,
    stringsAsFactors = FALSE
  )
}

result_df <- bind_rows(results)
print(result_df)
