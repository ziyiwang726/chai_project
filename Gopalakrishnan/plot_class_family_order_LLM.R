# Run install_packages.R once in a fresh R session before sourcing this script.
suppressPackageStartupMessages({
  library(ape)
  library(CATMicrobiome)
  library(coin)
  library(admix)
  library(mvtnorm)
  library(mixtools)
  library(phyloseq)
  library(mclust)
  library(adaptMT)
  library(AdaPTGMM)
  library(splines)
  library(IHW)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(tidyr)
  library(readr)
  library(auctestr)
  library(ggplot2)
  library(ggrepel)
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

get_worker_cores <- function() {
  env_cores <- suppressWarnings(as.integer(Sys.getenv("LLM_CORES", unset = "")))
  if (!is.na(env_cores) && env_cores > 0) {
    return(env_cores)
  }

  detected <- parallel::detectCores(logical = TRUE)
  if (is.na(detected) || detected < 1) {
    detected <- suppressWarnings(as.integer(
      tryCatch(system("getconf _NPROCESSORS_ONLN", intern = TRUE), error = function(e) NA_character_)
    ))
  }
  if (is.na(detected) || detected < 1) {
    detected <- 1L
  }

  max(1L, detected - 1L)
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

provider_dir_name <- function(provider) {
  switch(
    provider,
    openai = "openAIGenerated",
    claude = "claudeGenerated",
    gemini = "geminiGenerated",
    combined = "combined"
  )
}

provider_aux_filename <- function(provider, aux_source, aux_mode) {
  if (provider == "openai") {
    sprintf("d1Taxonomy_OTU_ABCD_%s_and_below_%s_openai.csv", aux_source, aux_mode)
  } else {
    sprintf("d1Taxonomy_OTU_ABCD_%s_and_below_%s.csv", aux_source, aux_mode)
  }
}

get_chai_lfdr <- function(chai_fit) {
  if (!is.null(chai_fit$lFDR)) return(chai_fit$lFDR)
  if (!is.null(chai_fit$clFDR)) return(chai_fit$clFDR)
  stop("CHAI result does not contain lFDR or clFDR.")
}

script_path <- get_script_path()
script_dir <- if (!is.na(script_path)) dirname(script_path) else normalizePath(getwd())
provider <- normalize_provider(Sys.getenv("LLM_PROVIDER", unset = "combined"))
enable_fdrreg <- tolower(Sys.getenv("LLM_ENABLE_FDRREG", unset = "false")) %in%
  c("1", "true", "yes")
worker_cores <- get_worker_cores()

shared_cache_root <- file.path(script_dir, "cache", "featureSelectionLines", "shared", "target_levels")
provider_cache_root <- file.path(script_dir, "cache", "featureSelectionLines", provider)
provider_plot_root <- file.path(script_dir, "plots", "featureSelectionLines", provider)
dir.create(shared_cache_root, recursive = TRUE, showWarnings = FALSE)
dir.create(provider_cache_root, recursive = TRUE, showWarnings = FALSE)
dir.create(provider_plot_root, recursive = TRUE, showWarnings = FALSE)

source(file.path(script_dir, "code", "chai.R"))
source(file.path(script_dir, "code", "color_helper.R"))
source(file.path(script_dir, "code", "conditionalParam.R"))
source(file.path(script_dir, "code", "naiveRemoveOneObs.R"))
source(file.path(script_dir, "code", "rGaussianMix.R"))
source(file.path(script_dir, "code", "utils.R"))
source(file.path(script_dir, "code", "performance.R"))

load_or_fit <- function(cache_file, fit_fun, valid_fun = function(x) TRUE) {
  if (file.exists(cache_file)) {
    obj <- readRDS(cache_file)
    if (isTRUE(valid_fun(obj))) return(obj)
  }
  obj <- fit_fun()
  saveRDS(obj, cache_file)
  obj
}

warm_shared_target_caches <- function(ps, target_ranks) {
  invisible(lapply(target_ranks, function(rank) {
    rank_norm <- tolower(rank)
    rank_ps <- load_or_fit(
      file.path(shared_cache_root, paste0("cache_", rank_norm, "_target_ps.rds")),
      function() build_rank_ps(ps, rank_norm)
    )

    cs <- colSums(otu_table(rank_ps))
    cs[cs == 0] <- 1
    t_abun <- t(sweep(otu_table(rank_ps), 2, cs, "/") * 100)
    outcome <- factor(sample_data(rank_ps)$BinOutcomes, levels = c("NR", "R"))

    load_or_fit(
      file.path(shared_cache_root, paste0("cache_", rank_norm, "_target_pz.rds")),
      function() compute_pz_wilcox(t_abun, outcome)
    )

    NULL
  }))
}

clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    str_replace("^.*?__", "") %>%
    str_replace_all("\\[|\\]", "") %>%
    str_squish()
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

build_ps <- function() {
  otu_path <- system.file("extdata", "d1OTUtable.csv", package = "CATMicrobiome")
  tax_path <- system.file("extdata", "d1Taxonomy.csv", package = "CATMicrobiome")
  meta_path <- system.file("extdata", "d1Meta.csv", package = "CATMicrobiome")
  tree_path <- system.file("extdata", "d1Tree.tree", package = "CATMicrobiome")

  otutable <- read.csv(otu_path, header = TRUE, row.names = 1)
  taxonomy <- read.csv(tax_path, header = TRUE, row.names = 1)
  metaData <- read.csv(meta_path, header = TRUE, row.names = 1)
  tree <- read.tree(tree_path)

  otus_keep <- intersect(tree$tip.label, rownames(taxonomy))
  tree <- keep.tip(tree, otus_keep)
  taxonomy <- taxonomy[match(tree$tip.label, rownames(taxonomy)), , drop = FALSE]
  otutable <- otutable[otus_keep, , drop = FALSE]

  phyloseq(
    otu_table(as.matrix(otutable), taxa_are_rows = TRUE),
    tax_table(as.matrix(taxonomy)),
    sample_data(metaData),
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

rank_display_name <- function(rank) {
  rank_norm <- tolower(rank)
  switch(rank_norm, otu = "OTU", tools::toTitleCase(rank_norm))
}

make_main_df <- function(t_abun, pz_df, ps_obj, rank) {
  if (tolower(rank) == "otu") {
    pz_map <- if ("feature" %in% names(pz_df)) {
      pz_df %>% transmute(
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
        mutate(rank_value = otu_rep, taxon_key = otu_rep)
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
    pz_df %>% transmute(
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

stouffer_nonzero <- function(v) {
  v2 <- v[!is.na(v) & v != 0]
  if (!length(v2)) return(0)
  stouffer_z(v2, ignore.na = FALSE)
}

make_aux_tax_level <- function(aux_otu, level, aux_covariates) {
  if (tolower(level) == "otu") {
    otu_col <- names(aux_otu)[tolower(names(aux_otu)) == "otuname"]
    if (length(otu_col) != 1) {
      stop("Could not uniquely find OTU column 'OTUname' in aux table.")
    }

    return(
      aux_otu %>%
        filter(!is.na(.data[[otu_col]]), .data[[otu_col]] != "") %>%
        group_by(.data[[otu_col]]) %>%
        summarise(
          across(all_of(aux_covariates), stouffer_nonzero, .names = "x{col}"),
          .groups = "drop"
        ) %>%
        transmute(
          otu_rep = .data[[otu_col]],
          !!!setNames(lapply(paste0("x", aux_covariates), function(col) rlang::expr(.data[[!!col]])),
                      paste0("x", aux_covariates)),
          taxon_key = otu_rep
        )
    )
  }

  level_col <- names(aux_otu)[tolower(names(aux_otu)) == tolower(level)]
  if (length(level_col) != 1) {
    stop("Could not uniquely find level column '", level, "' in aux table.")
  }

  aux_otu %>%
    filter(!is.na(.data[[level_col]]), .data[[level_col]] != "") %>%
    group_by(.data[[level_col]]) %>%
    summarise(
      across(all_of(aux_covariates), stouffer_nonzero, .names = "x{col}"),
      .groups = "drop"
    ) %>%
    mutate(taxon_key = clean_taxon(.data[[level_col]]))
}

resolve_aux_file <- function(provider, aux_source, aux_mode) {
  file.path(
    script_dir,
    "sideCov",
    provider_dir_name(provider),
    aux_mode,
    "selectionTargetOTU",
    provider_aux_filename(provider, aux_source, aux_mode)
  )
}

build_ns_formula <- function(X, requested_df) {
  covariate_unique_counts <- vapply(X, function(col) {
    length(unique(col[is.finite(col)]))
  }, integer(1))
  varying_covariates <- names(covariate_unique_counts)[covariate_unique_counts > 1]

  if (!length(varying_covariates)) {
    return(~ 1)
  }

  terms <- vapply(varying_covariates, function(var_name) {
    n_unique <- covariate_unique_counts[[var_name]]
    if (n_unique <= 2) {
      var_name
    } else {
      sprintf("ns(%s, df=%d)", var_name, min(requested_df, n_unique - 1))
    }
  }, character(1))

  as.formula(paste("~", paste(terms, collapse = " + ")))
}

compute_ihw_rejections <- function(p, covariate, covariate_name, q_levels) {
  finite_vals <- covariate[is.finite(covariate)]
  if (length(unique(finite_vals)) <= 1) {
    message("Skipping IHW for ", covariate_name, " because the covariate is constant.")
    return(rep(NA_real_, length(q_levels)))
  }

  tryCatch(
    sapply(q_levels, function(q) {
      ihw_fit <- ihw(pvalues = p, covariates = covariate, alpha = q, nbins = 5)
      rejections(ihw_fit)
    }),
    error = function(e) {
      warning("IHW failed for ", covariate_name, ". Error: ", conditionMessage(e))
      rep(NA_real_, length(q_levels))
    }
  )
}

fit_ordershape <- function(p, covariate, covariate_name) {
  if (!requireNamespace("OrderShapeEM", quietly = TRUE)) return(NULL)
  finite_vals <- covariate[is.finite(covariate)]
  if (length(unique(finite_vals)) <= 1) {
    message("Skipping OrderShapeEM for ", covariate_name, " because the covariate is constant.")
    return(NULL)
  }

  tryCatch(
    OrderShapeEM::OrderShapeEM(
      pvals = p,
      order.var = covariate,
      OrderShapeEM::OrderShapeEM.control(trace = TRUE)
    ),
    error = function(e) {
      warning("OrderShapeEM failed for ", covariate_name, ". Error: ", conditionMessage(e))
      NULL
    }
  )
}

plot_rejections_vs_q <- function(long_table, title) {
  ord <- method_order(long_table$Method)
  long_table <- long_table %>% mutate(Method = factor(Method, levels = ord))
  pal <- make_method_palette(levels(long_table$Method), K = 12)

  endpoints <- long_table %>%
    filter(is.finite(Rejections)) %>%
    group_by(Method) %>%
    filter(q == max(q, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    mutate(label = as.character(Method))

  ggplot(long_table, aes(x = q, y = Rejections, color = Method, group = Method)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 1.7, na.rm = TRUE) +
    ggrepel::geom_text_repel(
      data = endpoints,
      aes(label = label),
      nudge_x = 0.006,
      direction = "y",
      hjust = 0,
      segment.size = 0.2,
      box.padding = 0.15,
      point.padding = 0.1,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = pal, breaks = ord, drop = FALSE) +
    scale_x_continuous(
      breaks = sort(unique(long_table$q)),
      expand = expansion(mult = c(0.02, 0.18))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = title,
      x = "FDR threshold (q)",
      y = "Number of features selected",
      color = "Method"
    ) +
    theme_minimal() +
    theme(plot.margin = margin(5.5, 80, 5.5, 5.5))
}

run_combo <- function(ps, provider, rank, aux_source, aux_mode, q_levels) {
  message("Running provider: ", provider,
          " | mode: ", aux_mode,
          " | target: ", rank,
          " | sideCov: ", aux_source)

  rank_norm <- tolower(rank)
  aux_source_norm <- tolower(aux_source)
  mode_plot_root <- file.path(provider_plot_root, aux_mode)
  combo_cache_root <- file.path(provider_cache_root, aux_mode, paste0("aux_", aux_source_norm))
  dir.create(mode_plot_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(combo_cache_root, recursive = TRUE, showWarnings = FALSE)

  rank_ps <- load_or_fit(
    file.path(shared_cache_root, paste0("cache_", rank_norm, "_target_ps.rds")),
    function() build_rank_ps(ps, rank_norm)
  )

  cs <- colSums(otu_table(rank_ps))
  cs[cs == 0] <- 1
  t_abun <- t(sweep(otu_table(rank_ps), 2, cs, "/") * 100)
  outcome <- factor(sample_data(rank_ps)$BinOutcomes, levels = c("NR", "R"))

  pz <- load_or_fit(
    file.path(shared_cache_root, paste0("cache_", rank_norm, "_target_pz.rds")),
    function() compute_pz_wilcox(t_abun, outcome)
  )

  main_df <- make_main_df(t_abun = t_abun, pz_df = pz, ps_obj = rank_ps, rank = rank_norm)

  aux_file <- resolve_aux_file(provider, aux_source_norm, aux_mode)
  if (!file.exists(aux_file)) stop("Missing aux file: ", aux_file)
  aux_otu <- readr::read_csv(aux_file, show_col_types = FALSE)

  required_aux_covariates <- c("A", "B", "C")
  missing_aux_covariates <- setdiff(required_aux_covariates, names(aux_otu))
  if (length(missing_aux_covariates)) {
    stop("Auxiliary OTU CSV is missing required columns: ",
         paste(missing_aux_covariates, collapse = ", "))
  }
  aux_covariates <- intersect(c("A", "B", "C", "D"), names(aux_otu))
  model_covariates <- paste0("x", required_aux_covariates)

  aux_rank <- make_aux_tax_level(aux_otu, level = rank_norm, aux_covariates = aux_covariates)

  df_model <- load_or_fit(
    file.path(combo_cache_root, paste0("cache_df_model_target_", rank_norm, "_aux_", aux_source_norm, ".rds")),
    function() {
      main_df %>%
        left_join(aux_rank %>% select(taxon_key, all_of(model_covariates)), by = "taxon_key") %>%
        mutate(across(all_of(model_covariates), ~ tidyr::replace_na(.x, 0)))
    }
  )

  z <- df_model$z_main
  p <- df_model$p_main
  X <- as.data.frame(df_model[, model_covariates, drop = FALSE])
  keep <- complete.cases(z, p, X)
  z <- z[keep]
  p <- p[keep]
  X <- X[keep, , drop = FALSE]

  if (length(z) < 10) stop("Too few complete rows for target rank ", rank, ".")

  p_adjusted <- p.adjust(p, method = "BH")
  bh_sel <- which(p_adjusted <= 0.05)

  chai_fit <- tryCatch(
    load_or_fit(
      file.path(combo_cache_root, paste0("cache_chai_target_", rank_norm, "_aux_", aux_source_norm, ".rds")),
      function() {
        set.seed(123)
        chai(z, X, B = 100)
      },
      valid_fun = function(obj) {
        lfdr <- tryCatch(get_chai_lfdr(obj), error = function(e) NULL)
        !is.null(lfdr) && length(lfdr) == length(z)
      }
    ),
    error = function(e) {
      warning("CHAI failed for target ", rank, " / aux ", aux_source, ". Error: ",
              conditionMessage(e))
      NULL
    }
  )
  chai_lfdr <- if (!is.null(chai_fit)) get_chai_lfdr(chai_fit) else rep(NA_real_, length(z))

  formula_ns <- lapply(c(2, 4, 6), function(df_value) build_ns_formula(X, df_value))

  glm_ns <- tryCatch(
    load_or_fit(
      file.path(combo_cache_root, paste0("cache_glm_target_", rank_norm, "_aux_", aux_source_norm, ".rds")),
      function() adapt_glm(X, pvals = p, alphas = q_levels,
                           pi_formulas = formula_ns, mu_formulas = formula_ns),
      valid_fun = function(obj) length(obj$nrejs) >= length(q_levels)
    ),
    error = function(e) {
      warning("adapt_glm failed for target ", rank, " / aux ", aux_source, ". Error: ",
              conditionMessage(e))
      NULL
    }
  )

  gmm_ns_p <- tryCatch(
    load_or_fit(
      file.path(combo_cache_root, paste0("cache_gmm_p_target_", rank_norm, "_aux_", aux_source_norm, ".rds")),
      function() adapt_gmm(X, pvals = p, alphas = q_levels, beta_formulas = formula_ns),
      valid_fun = function(obj) length(obj$nrejs) >= length(q_levels)
    ),
    error = function(e) {
      warning("adapt_gmm(p) failed for target ", rank, " / aux ", aux_source, ". Error: ",
              conditionMessage(e))
      NULL
    }
  )

  gmm_ns_z <- tryCatch(
    load_or_fit(
      file.path(combo_cache_root, paste0("cache_gmm_z_target_", rank_norm, "_aux_", aux_source_norm, ".rds")),
      function() adapt_gmm(X, z = z, alphas = q_levels,
                           beta_formulas = formula_ns, testing = "two_sided"),
      valid_fun = function(obj) length(obj$nrejs) >= length(q_levels)
    ),
    error = function(e) {
      warning("adapt_gmm(z) failed for target ", rank, " / aux ", aux_source, ". Error: ",
              conditionMessage(e))
      NULL
    }
  )

  os <- tryCatch(
    load_or_fit(
      file.path(combo_cache_root, paste0("cache_ordershape_target_", rank_norm, "_aux_", aux_source_norm, ".rds")),
      function() {
        list(
          x1 = fit_ordershape(p, X[, 1], colnames(X)[1]),
          x2 = fit_ordershape(p, X[, 2], colnames(X)[2]),
          x3 = fit_ordershape(p, X[, 3], colnames(X)[3])
        )
      }
    ),
    error = function(e) {
      warning("OrderShapeEM failed for target ", rank, " / aux ", aux_source, ". Error: ",
              conditionMessage(e))
      NULL
    }
  )

  fdr_theo <- if (enable_fdrreg && requireNamespace("FDRreg", quietly = TRUE)) {
    tryCatch(
      load_or_fit(
        file.path(combo_cache_root, paste0("cache_fdrreg_target_", rank_norm, "_aux_", aux_source_norm, ".rds")),
        function() {
          set.seed(123)
          FDRreg::FDRreg(z, as.matrix(X), nulltype = "theoretical")
        }
      ),
      error = function(e) {
        warning("FDRreg failed for target ", rank, " / aux ", aux_source, ". Error: ",
                conditionMessage(e))
        NULL
      }
    )
  } else {
    NULL
  }

  chai_rejs <- if (!is.null(chai_fit)) {
    sapply(q_levels, function(q) length(lFDRselect(chai_lfdr, q, 1)))
  } else {
    rep(NA_real_, length(q_levels))
  }
  adapt_glm_ns_rej <- if (!is.null(glm_ns)) glm_ns$nrejs[seq_along(q_levels)] else rep(NA_real_, length(q_levels))
  adapt_gmm_p_ns_rej <- if (!is.null(gmm_ns_p)) gmm_ns_p$nrejs[seq_along(q_levels)] else rep(NA_real_, length(q_levels))
  adapt_gmm_z_ns_rej <- if (!is.null(gmm_ns_z)) gmm_ns_z$nrejs[seq_along(q_levels)] else rep(NA_real_, length(q_levels))

  ordershapeem_x1_rejs <- if (!is.null(os) && !is.null(os$x1)) {
    sapply(q_levels, function(q) sum(os$x1$fdr <= q))
  } else {
    rep(NA_real_, length(q_levels))
  }
  ordershapeem_x2_rejs <- if (!is.null(os) && !is.null(os$x2)) {
    sapply(q_levels, function(q) sum(os$x2$fdr <= q))
  } else {
    rep(NA_real_, length(q_levels))
  }
  ordershapeem_x3_rejs <- if (!is.null(os) && !is.null(os$x3)) {
    sapply(q_levels, function(q) sum(os$x3$fdr <= q))
  } else {
    rep(NA_real_, length(q_levels))
  }

  rdrreg_rejs <- if (!is.null(fdr_theo)) {
    sapply(q_levels, function(q) length(which(fdr_theo$FDR <= q)))
  } else {
    rep(NA_real_, length(q_levels))
  }

  ihw_1_rejs <- compute_ihw_rejections(p, X[, 1], colnames(X)[1], q_levels)
  ihw_2_rejs <- compute_ihw_rejections(p, X[, 2], colnames(X)[2], q_levels)
  ihw_3_rejs <- compute_ihw_rejections(p, X[, 3], colnames(X)[3], q_levels)
  bh_rejs <- sapply(q_levels, function(q) length(which(p_adjusted <= q)))

  rank_df <- data.frame(
    q = q_levels,
    chai = chai_rejs,
    adapt_glm = adapt_glm_ns_rej,
    adapt_gmm_p = adapt_gmm_p_ns_rej,
    adapt_gmm_z = adapt_gmm_z_ns_rej,
    OrderShapeEM_x1 = ordershapeem_x1_rejs,
    OrderShapeEM_x2 = ordershapeem_x2_rejs,
    OrderShapeEM_x3 = ordershapeem_x3_rejs,
    FDRreg = rdrreg_rejs,
    IHW_x1 = ihw_1_rejs,
    IHW_x2 = ihw_2_rejs,
    IHW_x3 = ihw_3_rejs,
    BH = bh_rejs
  )

  rank_long <- rank_df %>%
    pivot_longer(cols = -q, names_to = "Method", values_to = "Rejections") %>%
    mutate(Method = as.character(Method))

  out_pdf <- file.path(
    mode_plot_root,
    sprintf("feature_selection_target_%s_aux_%s_%s.pdf", rank_norm, aux_source_norm, aux_mode)
  )
  p_rank <- plot_rejections_vs_q(
    rank_long,
    title = paste0(
      "Feature Selection vs FDR (",
      provider,
      ", ",
      aux_mode,
      ", target=",
      rank_display_name(rank),
      ", sideCov=",
      rank_display_name(aux_source_norm),
      ")"
    )
  )
  ggsave(out_pdf, p_rank, width = 10, height = 5)

  tibble::tibble(
    provider = provider,
    aux_mode = aux_mode,
    aux_source = aux_source_norm,
    target_rank = rank_norm,
    aux_file = normalizePath(aux_file, mustWork = FALSE),
    output_pdf = normalizePath(out_pdf, mustWork = FALSE),
    gaussian_clusters = if (!is.null(chai_fit) && !is.null(chai_fit$K)) as.integer(chai_fit$K[[1]]) else NA_integer_,
    n_features = length(z),
    chai_selected_q05 = if (!is.null(chai_fit)) length(lFDRselect(chai_lfdr, 0.05, 1)) else NA_integer_,
    chai_selected_q10 = if (!is.null(chai_fit)) length(lFDRselect(chai_lfdr, 0.10, 1)) else NA_integer_,
    chai_selected_q15 = if (!is.null(chai_fit)) length(lFDRselect(chai_lfdr, 0.15, 1)) else NA_integer_,
    chai_selected_q20 = if (!is.null(chai_fit)) length(lFDRselect(chai_lfdr, 0.20, 1)) else NA_integer_,
    bh_selected_q05 = length(bh_sel)
  )
}

ps <- build_ps()
q_levels <- seq(0.01, 0.20, by = 0.01)
default_aux_modes <- c("unweighted", "weighted")
default_aux_sources <- c("class", "order", "family", "genus", "species")
default_target_ranks <- c("family", "genus", "species", "otu")

aux_modes <- parse_env_list("LLM_AUX_MODES", default_aux_modes)
aux_sources <- parse_env_list("LLM_AUX_SOURCES", default_aux_sources)
target_ranks <- parse_env_list("LLM_TARGET_RANKS", default_target_ranks)
warm_shared_target_caches(ps, target_ranks)

message("Using ", worker_cores, " worker process(es) for feature-selection line plots.")
run_grid <- expand.grid(
  aux_mode = aux_modes,
  aux_source = aux_sources,
  rank = target_ranks,
  stringsAsFactors = FALSE
)

run_one_combo <- function(i) {
  row <- run_grid[i, ]
  run_combo(
    ps = ps,
    provider = provider,
    rank = row$rank,
    aux_source = row$aux_source,
    aux_mode = row$aux_mode,
    q_levels = q_levels
  )
}

if (.Platform$OS.type == "unix" && worker_cores > 1) {
  summary_rows <- parallel::mclapply(
    seq_len(nrow(run_grid)),
    run_one_combo,
    mc.cores = worker_cores,
    mc.preschedule = FALSE
  )
} else {
  summary_rows <- lapply(seq_len(nrow(run_grid)), run_one_combo)
}

summary_df <- bind_rows(summary_rows)
summary_file <- file.path(provider_plot_root, "chai_cluster_summary.csv")
if (file.exists(summary_file)) {
  existing_summary <- readr::read_csv(summary_file, show_col_types = FALSE)
  summary_df <- bind_rows(existing_summary, summary_df) %>%
    distinct(provider, aux_mode, aux_source, target_rank, .keep_all = TRUE)
}
summary_df <- summary_df %>%
  arrange(aux_mode, target_rank, aux_source)
readr::write_csv(summary_df, summary_file)
print(summary_df)
message("Saved CHAI cluster summary: ", summary_file)
