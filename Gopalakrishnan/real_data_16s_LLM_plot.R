# Run install_packages.R once in a fresh R session before sourcing this script.
suppressPackageStartupMessages({
  library(devtools)
  library(ape)
  library(CATMicrobiome)
  library(coin)
  library(locfdr)
  library(admix)
  library(mvtnorm)
  library(mixtools)
  library(phyloseq)
  library(KScorrect)
  library(mclust)
  library(adaptMT)
  library(AdaPTGMM)
  library(splines2)
  library(splines)
  library(IHW)
  library(FDRreg)
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

provider_aux_filename <- function(provider, side_mode) {
  switch(
    provider,
    openai = sprintf("d1Taxonomy_OTU_ABCD_genus_and_below_%s_openai.csv", side_mode),
    claude = sprintf("d1Taxonomy_OTU_ABCD_genus_and_below_%s.csv", side_mode),
    gemini = sprintf("d1Taxonomy_OTU_ABCD_genus_and_below_%s.csv", side_mode),
    combined = sprintf("d1Taxonomy_OTU_ABCD_genus_and_below_%s.csv", side_mode)
  )
}

get_chai_lfdr <- function(chai_fit) {
  if (!is.null(chai_fit$lFDR)) return(chai_fit$lFDR)
  if (!is.null(chai_fit$clFDR)) return(chai_fit$clFDR)
  stop("CHAI result does not contain lFDR or clFDR.")
}

script_path <- get_script_path()
script_dir <- if (!is.na(script_path)) dirname(script_path) else normalizePath(getwd())
provider <- normalize_provider(Sys.getenv("LLM_PROVIDER", unset = "combined"))
side_mode <- match.arg(tolower(Sys.getenv("LLM_SIDE_MODE", unset = "unweighted")),
                       c("unweighted", "weighted"))
enable_fdrreg <- tolower(Sys.getenv("LLM_ENABLE_FDRREG", unset = "true")) %in%
  c("1", "true", "yes")

default_aux_csv <- file.path(
  script_dir,
  "sideCov",
  provider_dir_name(provider),
  side_mode,
  "selectionTargetOTU",
  provider_aux_filename(provider, side_mode)
)
legacy_aux_csv <- file.path(
  script_dir,
  "sideCov",
  "Gemini3Generated",
  "unweighted",
  "d1Taxonomy_OTU_ABCD_genus_and_below.csv"
)
aux_csv <- Sys.getenv("LLM_AUX_OTU_CSV", unset = "")
if (!nzchar(aux_csv)) {
  aux_csv <- if (file.exists(default_aux_csv)) default_aux_csv else legacy_aux_csv
}
if (!file.exists(aux_csv)) {
  stop("Auxiliary OTU CSV not found: ", aux_csv)
}

cache_root <- Sys.getenv(
  "LLM_CACHE_ROOT",
  unset = file.path(script_dir, "cache", "real_data_16s_LLM_plotYS", provider, side_mode)
)
plot_root <- Sys.getenv(
  "LLM_PLOT_ROOT",
  unset = file.path(script_dir, "plots", "real_data_16s_LLM_plotYS", provider, side_mode)
)
dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)
cache_file <- function(name) file.path(cache_root, name)
plot_file <- function(name) file.path(plot_root, name)
save_plot <- function(name, plot_obj, width = 10, height = 5) {
  ggsave(plot_file(name), plot_obj, width = width, height = height)
}
run_label <- sprintf("%s (%s)", provider, side_mode)
message("Running real_data_16s_LLM_plotYS.R with provider = ", provider,
        ", side mode = ", side_mode,
        ", aux CSV = ", normalizePath(aux_csv, mustWork = FALSE))

source(file.path(script_dir, "code", "chai.R"))
source(file.path(script_dir, "code", "color_helper.R"))
source(file.path(script_dir, "code", "conditionalParam.R"))
source(file.path(script_dir, "code", "naiveRemoveOneObs.R"))
source(file.path(script_dir, "code", "rGaussianMix.R"))
source(file.path(script_dir, "code", "utils.R"))
source(file.path(script_dir, "code", "performance.R"))
otuPath <- system.file("extdata","d1OTUtable.csv",
                       package = "CATMicrobiome")
otutable <- read.csv(otuPath,header=TRUE,row.names = 1)

taxonomyPathLocal <- file.path(script_dir, "LLMCode", "d1Taxonomy.csv")
taxonomyPath <- if (file.exists(taxonomyPathLocal)) {
  taxonomyPathLocal
} else {
  system.file("extdata","d1Taxonomy.csv", package = "CATMicrobiome")
}
taxonomy <- read.csv(taxonomyPath,header=TRUE,row.names = 1)

metaPath <- system.file("extdata","d1Meta.csv",
                        package = "CATMicrobiome")
metaData <- read.csv(metaPath,header=TRUE,row.names = 1)

treePath <- system.file("extdata","d1Tree.tree",
                        package = "CATMicrobiome")
tree <- read.tree(treePath)

otus_tree   <- tree$tip.label       # 1437
otus_tax    <- rownames(taxonomy)   # 1455
otus_keep   <- intersect(otus_tree, otus_tax)   # Keep 1437

# Prune tree and taxonomy
tree <- keep.tip(tree, otus_keep)
taxonomy <- taxonomy[otus_keep, , drop = FALSE]
taxonomy <- taxonomy[match(tree$tip.label, rownames(taxonomy)), , drop = FALSE]  # reorder to tree tips
otutable <- otutable[otus_keep, , drop = FALSE]
# 1437

OTU <- otu_table(as.matrix(otutable), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(taxonomy))
SAM <- sample_data(metaData)
TR  <- phy_tree(tree)

ps <- phyloseq(OTU, TAX, SAM, TR)

###################### Genus level ps #######################
if (file.exists(cache_file("cache_genus_ps.rds"))) {
  genus <- readRDS(cache_file("cache_genus_ps.rds"))
} else {
  genus <- tax_glom(ps, taxrank='genus')
  saveRDS(genus, cache_file("cache_genus_ps.rds"))
}


###################### Abundance table #######################
# genus
cs <- colSums(otu_table(genus))
cs[cs == 0] <- 1
genus_abun <- sweep(otu_table(genus), 2, cs, "/") * 100
t_genus_abun <- t(genus_abun)

########################## Outcome ###################################
# Check the metadata order matching to the abundance subject order
rownames(sample_data(genus)) == rownames(t_genus_abun)

outcome <- factor(sample_data(genus)$BinOutcomes, levels = c("NR","R"))
table(outcome)


##################### p-value/z-value ###########################
compute_pz_wilcox <- function(t_abun, outcome) {
  stopifnot(nrow(t_abun) == length(outcome))

  p_z_values <- sapply(colnames(t_abun), function(feature) {
    features <- t_abun[, feature, drop = FALSE]
    wilcox <- wilcox_test(features[, 1] ~ outcome)

    c(
      p_value = as.numeric(pvalue(wilcox)),
      z_value = as.numeric(statistic(wilcox))
    )
  }, simplify = "matrix")
  p_z_values <- as.data.frame(t(p_z_values))
  p_z_values$feature <- rownames(p_z_values)
  rownames(p_z_values) <- NULL

  return(p_z_values)
}

# genus's p-value and z-value
if (file.exists(cache_file("cache_genus_pz.rds"))) {
  genus_pz <- readRDS(cache_file("cache_genus_pz.rds"))
} else {
  genus_pz <- compute_pz_wilcox(t_genus_abun, outcome)
  saveRDS(genus_pz, cache_file("cache_genus_pz.rds"))
}


#################### Build OTU for each level #####################
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)



clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    str_replace("^.*?__", "") %>%        # Remove k__/f__/g__
    str_replace_all("\\[|\\]", "") %>%   # Remove []
    str_squish()                         # Remove extra space
}

make_main_df <- function(t_abun, pz_df, ps_obj, rank) {
  # rank: the tax_table column you want to join on, e.g. "family" / "genus" / "order" / "class"

  tax_df <- as.data.frame(tax_table(ps_obj)) %>%
    rownames_to_column("otu_rep")

  # find the taxonomy column robustly (case-insensitive)
  rank_col <- names(tax_df)[tolower(names(tax_df)) == tolower(rank)]
  if (length(rank_col) != 1) {
    stop("Could not uniquely find taxonomy column for rank = '", rank, "'. ",
         "Available columns: ", paste(names(tax_df), collapse = ", "))
  }

  tax_map <- tax_df %>%
    transmute(
      otu_rep,
      rank_value = as.character(.data[[rank_col]])
    )

  # p/z mapping (supports either pz_df$feature or assumes same order as colnames(t_abun))
  pz_map <- if ("feature" %in% names(pz_df)) {
    pz_df %>%
      transmute(
        otu_rep = feature,
        p_main  = as.numeric(p_value),
        z_main  = as.numeric(z_value)
      )
  } else {
    data.frame(
      otu_rep = colnames(t_abun),
      p_main  = as.numeric(pz_df$p_value),
      z_main  = as.numeric(pz_df$z_value),
      stringsAsFactors = FALSE
    )
  }

  main_df <- pz_map %>%
    left_join(tax_map, by = "otu_rep") %>%
    mutate(taxon_key = clean_taxon(rank_value))

  main_df
}

genus_pz_tax = make_main_df(t_abun = t_genus_abun,  pz_df = genus_pz,
                            ps_obj = genus,  rank = "genus")


#################### LLM results at family level ######################
library(dplyr)
library(readr)
library(auctestr)

genus_aux_otu <- read_csv(aux_csv, show_col_types = FALSE)
required_aux_covariates <- c("A", "B", "C")
missing_aux_covariates <- setdiff(required_aux_covariates, names(genus_aux_otu))
if (length(missing_aux_covariates)) {
  stop("Auxiliary OTU CSV is missing required columns: ",
       paste(missing_aux_covariates, collapse = ", "))
}
aux_covariates <- intersect(c("A", "B", "C", "D"), names(genus_aux_otu))
model_covariates <- paste0("x", required_aux_covariates)

# Aggregation using stouffer_z, only stouffer non-zero values
stouffer_nonzero <- function(v) {
  v2 <- v[!is.na(v) & v != 0]
  if (length(v2) == 0) return(0)
  stouffer_z(v2, ignore.na = FALSE)
}

make_aux_tax_level <- function(aux_otu,
                               level,
                               agg_z_method = c("stouffer", "median")) {
  agg_z_method <- match.arg(agg_z_method)

  # pick aggregation function
  agg_fun <- switch(
    agg_z_method,
    stouffer = stouffer_nonzero,
    median   = median_nonzero
  )

  aux_lvl <- aux_otu %>%
    filter(!is.na(.data[[level]]), .data[[level]] != "") %>%
    group_by(.data[[level]]) %>%
    summarise(
      across(
        .cols  = all_of(aux_covariates),
        .fns   = agg_fun,
        .names = "x{col}"
      ),
      .groups = "drop"
    ) %>%
    mutate(taxon_key = clean_taxon(.data[[level]]))

  aux_lvl
}

# Genus - Used genus and below, aggregated at genus level
gen_ag_gen <- make_aux_tax_level(genus_aux_otu, level = "genus",
                                 agg_z_method = "stouffer")



############################ Fit Models ###############################
if (file.exists(cache_file("cache_df_model.rds"))) {
  df_model <- readRDS(cache_file("cache_df_model.rds"))
} else {
  df_model <- genus_pz_tax %>%
    dplyr::left_join(
      gen_ag_gen %>% dplyr::select(taxon_key, dplyr::all_of(model_covariates)),
      by = "taxon_key"
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(model_covariates), ~ tidyr::replace_na(.x, 0)))
  saveRDS(df_model, cache_file("cache_df_model.rds"))
}

z <- df_model$z_main
p <- df_model$p_main
X <- as.data.frame(df_model[, model_covariates, drop = FALSE])

############################ Visualization ###############################
df <- data.frame(z, X)
df_long <- pivot_longer(df, cols = starts_with("x"), names_to = "variable", values_to = "value")
p_relationship <- ggplot(df_long, aes(x = value, y = z)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~variable, scales = "free_x") +
  theme_minimal() +
  labs(title = paste("Relationship between z and predictors", run_label))
save_plot("16s_wilcoxz_LLM_z_vs_predictors.pdf", p_relationship, width = 10, height = 5)

library(GGally)
p_llm_z <- ggpairs(df)
save_plot("16s_wilcoxz_LLM_pairs.pdf", p_llm_z, width = 10, height = 10)

################## BH #######################
p_adjusted <- p.adjust(p, method = "BH")
bh_sel <- which(p_adjusted <= 0.05)

################## chai #######################
if (file.exists(cache_file("cache_chai_genus_llm.rds"))) {
  chai_genus_llm <- readRDS(cache_file("cache_chai_genus_llm.rds"))
} else {
  set.seed(123)
  chai_genus_llm <- chai(z, X, B = 100)
  saveRDS(chai_genus_llm, cache_file("cache_chai_genus_llm.rds"))
}
chai_lfdr <- get_chai_lfdr(chai_genus_llm)
chai_clusters <- if (!is.null(chai_genus_llm$K)) as.integer(chai_genus_llm$K[[1]]) else NA_integer_
chai_sel_q05 <- lFDRselect(chai_lfdr, 0.05, 1)
chai_sel_q10 <- lFDRselect(chai_lfdr, 0.1, 1)
length(chai_sel_q05)
length(chai_sel_q10)

################## Adapt #######################
alphas <- seq(0.01, 0.20, by = 0.01)

# GLM
# ns
x <- data.frame(X)
covariate_unique_counts <- vapply(x, function(col) {
  length(unique(col[is.finite(col)]))
}, integer(1))
varying_covariates <- names(covariate_unique_counts)[covariate_unique_counts > 1]

build_ns_formula <- function(requested_df) {
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

formula_ns <- lapply(c(2, 4, 6), build_ns_formula)

if (file.exists(cache_file("cache_glm_ns_q20.rds"))) {
  glm_ns <- readRDS(cache_file("cache_glm_ns_q20.rds"))
} else {
  glm_ns <- adapt_glm(x, pvals = p, alphas=alphas,
                      pi_formulas = formula_ns, mu_formulas = formula_ns)
  saveRDS(glm_ns, cache_file("cache_glm_ns_q20.rds"))
}

# GMM - p
# ns
if (file.exists(cache_file("cache_gmm_ns_p_q20.rds"))) {
  gmm_ns_p <- readRDS(cache_file("cache_gmm_ns_p_q20.rds"))
} else {
  gmm_ns_p <- adapt_gmm(x, pvals = p,  alphas=alphas,
                        beta_formulas = formula_ns)
  saveRDS(gmm_ns_p, cache_file("cache_gmm_ns_p_q20.rds"))
}

# GMM - Z
if (file.exists(cache_file("cache_gmm_ns_z_q20.rds"))) {
  gmm_ns_z <- readRDS(cache_file("cache_gmm_ns_z_q20.rds"))
} else {
  gmm_ns_z <- adapt_gmm(x, z = z, alphas=alphas,
                        beta_formulas = formula_ns, testing = "two_sided")
  saveRDS(gmm_ns_z, cache_file("cache_gmm_ns_z_q20.rds"))
}
# alpha = 0.1: FDPhat 0.075, Number of Rej. 0
# alpha = 0.05: FDPhat 0.05, Number of Rej. 0


################## OrderShapeEM #######################
require(OrderShapeEM)
if (file.exists(cache_file("cache_ordershape.rds"))) {
  os_cache <- readRDS(cache_file("cache_ordershape.rds"))
  ordershape_x1 <- os_cache$x1
  ordershape_x2 <- os_cache$x2
  ordershape_x3 <- os_cache$x3
} else {
  ordershape_x1 <- OrderShapeEM(pvals = p, order.var = X[,1],
                                OrderShapeEM.control(trace = TRUE))
  ordershape_x2 <- OrderShapeEM(pvals = p, order.var = X[,2],
                                OrderShapeEM.control(trace = TRUE))
  ordershape_x3 <- OrderShapeEM(pvals = p, order.var = X[,3],
                                OrderShapeEM.control(trace = TRUE))
  saveRDS(list(x1 = ordershape_x1, x2 = ordershape_x2, x3 = ordershape_x3),
          cache_file("cache_ordershape.rds"))
}
sum(ordershape_x1$fdr <= 0.05)
sum(ordershape_x2$fdr <= 0.05)
sum(ordershape_x3$fdr <= 0.05)


######################## FDRreg ############################
fdreg_available <- enable_fdrreg && requireNamespace("FDRreg", quietly = TRUE)
if (fdreg_available) {
  if (file.exists(cache_file("cache_fdr_theo.rds"))) {
    fdr_theo <- readRDS(cache_file("cache_fdr_theo.rds"))
  } else {
    set.seed(123)
    fdr_theo <- FDRreg::FDRreg(z, as.matrix(X), nulltype = 'theoretical')
    saveRDS(fdr_theo, cache_file("cache_fdr_theo.rds"))
  }
  length(which(fdr_theo$FDR <= 0.05))
} else {
  message("FDRreg disabled or not available; skipping.")
  fdr_theo <- NULL
}

######################## IHW ############################
set.seed(123)
# ihw_x1_01 <- ihw(pvalues = p, covariates = X[,1], alpha = 0.01, nbins = 5)
# ihw_x1_02 <- ihw(pvalues = p, covariates = X[,1], alpha = 0.02, nbins = 5)
# ihw_x1_03 <- ihw(pvalues = p, covariates = X[,1], alpha = 0.03, nbins = 5)
# ihw_x1_04 <- ihw(pvalues = p, covariates = X[,1], alpha = 0.04, nbins = 5)
# ihw_x1_05 <- ihw(pvalues = p, covariates = X[,1], alpha = 0.05, nbins = 5)
# ihw_x1_06 <- ihw(pvalues = p, covariates = X[,1], alpha = 0.06, nbins = 5)
# ihw_x1_07 <- ihw(pvalues = p, covariates = X[,1], alpha = 0.07, nbins = 5)
# ihw_x1_08 <- ihw(pvalues = p, covariates = X[,1], alpha = 0.08, nbins = 5)
# ihw_x1_09 <- ihw(pvalues = p, covariates = X[,1], alpha = 0.09, nbins = 5)
# ihw_x1_1 <- ihw(pvalues = p, covariates = X[,1], alpha = 0.1, nbins = 5)
#
# ihw_x2_01 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.01, nbins = 5)
# ihw_x2_02 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.02, nbins = 5)
# ihw_x2_03 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.03, nbins = 5)
# ihw_x2_04 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.04, nbins = 5)
# ihw_x2_05 <- ihw(pvalues = p, covariates = X[,2], alpha = 0.05, nbins = 5)
# ihw_x2_06 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.06, nbins = 5)
# ihw_x2_07 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.07, nbins = 5)
# ihw_x2_08 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.08, nbins = 5)
# ihw_x2_09 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.09, nbins = 5)
# ihw_x2_1 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.1, nbins = 5)
#
# rejections(ihw_x1_05)  # 1   1
# rejections(ihw_x1_06)  # 1   1
# rejections(ihw_x1_07)  # 1   1
# rejections(ihw_x1_08)  # 0   1
# rejections(ihw_x1_09)  # 0   1
# rejections(ihw_x1_1)   # 0   2

###########################################################
# Summarize into a different q level table
q_levels <- alphas

# chai
chai_rejs <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_lfdr, q, 1))
})

# adapt_glm
adapt_glm_ns_rej <- glm_ns$nrejs[seq_along(q_levels)]

# adapt_gmm_p
adapt_gmm_p_ns_rej <- gmm_ns_p$nrejs[seq_along(q_levels)]
# adapt_gmm_p_bs_rej <- gmm_bs_p$nrejs[1:10]

# adapt_gmm_z
adapt_gmm_z_ns_rej <- gmm_ns_z$nrejs[seq_along(q_levels)]

# OrderShapeEM
ordershapeem_x1_rejs <- sapply(q_levels, function(q) {
  sum(ordershape_x1$fdr <= q)
})

ordershapeem_x2_rejs <- sapply(q_levels, function(q) {
  sum(ordershape_x2$fdr <= q)
})

ordershapeem_x3_rejs <- sapply(q_levels, function(q) {
  sum(ordershape_x3$fdr <= q)
})

# FDRreg
rdrreg_rejs <- if (fdreg_available) {
  sapply(q_levels, function(q) length(which(fdr_theo$FDR <= q)))
} else {
  rep(NA, length(q_levels))
}

# IHW
# IHW
compute_ihw_rejections <- function(covariate, covariate_name) {
  finite_vals <- covariate[is.finite(covariate)]
  if (length(unique(finite_vals)) <= 1) {
    message("Skipping IHW for ", covariate_name, " because the covariate is constant.")
    return(rep(NA_integer_, length(q_levels)))
  }

  sapply(q_levels, function(q) {
    ihw_fit <- ihw(pvalues = p, covariates = covariate, alpha = q, nbins = 5)
    rejections(ihw_fit)
  })
}

ihw_1_rejs <- compute_ihw_rejections(X[,1], colnames(X)[1])
ihw_2_rejs <- compute_ihw_rejections(X[,2], colnames(X)[2])
ihw_3_rejs <- compute_ihw_rejections(X[,3], colnames(X)[3])

# BH
p_adjusted <- p.adjust(p, method = "BH")

bh_rejs <- sapply(q_levels, function(q) {
  length(which(p_adjusted <= q))
})

# combine into data frame
LLM_gen_genandbelow <- data.frame(
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

chai_summary <- tibble::tibble(
  provider = provider,
  side_mode = side_mode,
  aux_csv = normalizePath(aux_csv, mustWork = FALSE),
  gaussian_clusters = chai_clusters,
  n_features = length(z),
  chai_selected_q05 = length(chai_sel_q05),
  chai_selected_q10 = length(chai_sel_q10),
  bh_selected_q05 = length(bh_sel)
)
readr::write_csv(chai_summary, plot_file("chai_run_summary.csv"))
readr::write_csv(LLM_gen_genandbelow, plot_file("rejections_by_q.csv"))




###########################################################
# Plots
library(ggplot2)
library(tidyr)
library(ggrepel)
library(viridisLite)

# into long format
LLM_gen_genandbelow_long <- LLM_gen_genandbelow %>%
  pivot_longer(
    cols = -q,
    names_to = "Method",
    values_to = "Rejections"
  ) %>%
  mutate(Method = as.character(Method))

# order <- c("OrderShapeEM_x1", "OrderShapeEM_x2", "OrderShapeEM_x3",
#            "OrderShapeEM_x4", "OrderShapeEM_x5", "IHW_x1", "IHW_x2",
#            "IHW_x3", "IHW_x4", "IHW_x5",
#            "adapt_glm", "adapt_gmm_z", "adapt_gmm_p",
#            "FDRreg", "BH", "chai")


plot_rejections_vs_q <- function(long_table, title = "Number of Features Selected vs FDR Threshold",
                                 K_variants = 12, x_breaks = NULL) {

  ord <- method_order(long_table$Method)
  long_table <- long_table %>% mutate(Method = factor(Method, levels = ord))
  pal <- make_method_palette(levels(long_table$Method), K = K_variants)

  endpoints <- long_table %>%
    group_by(Method) %>%
    filter(q == max(q, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    mutate(label = as.character(Method))

  x_rng   <- range(long_table$q, na.rm = TRUE)
  x_nudge <- diff(x_rng) * 0.03

  p <- ggplot(long_table, aes(x = q, y = Rejections, group = Method, color = Method)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 1.8, na.rm = TRUE) +
    ggrepel::geom_text_repel(
      data = endpoints,
      aes(label = label),
      nudge_x = x_nudge,
      hjust = 0,
      direction = "y",
      segment.size = 0.2,
      box.padding = 0.15,
      point.padding = 0.10,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = pal, breaks = ord) +
    labs(title = title, x = "FDR threshold (q)", y = "Number of features selected") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.margin = margin(5.5, 60, 5.5, 5.5)
    )

  if (!is.null(x_breaks)) {
    p <- p + scale_x_continuous(breaks = x_breaks, expand = expansion(mult = c(0.02, 0.18)))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0.02, 0.18)))
  }

  p
}

p_LLM <- plot_rejections_vs_q(LLM_gen_genandbelow_long,
                                  title = paste("Number of Features Selected vs FDR Threshold (LLM, genus)", run_label),
                                  x_breaks = sort(unique(LLM_gen_genandbelow_long$q)))

save_plot("16s_wilcoxz_LLM_gen_genandbelow.pdf", p_LLM, width = 10, height = 5)



# Visualization related to chai results
chai_sel <- chai_sel_q05
df <- data.frame(z, X) %>%
  dplyr::mutate(
    grp = dplyr::case_when(
      dplyr::row_number() %in% bh_sel & dplyr::row_number() %in% chai_sel ~ "both",
      dplyr::row_number() %in% bh_sel ~ "BH_only",
      dplyr::row_number() %in% chai_sel ~ "chai_only",
      TRUE ~ "neither"
    )
  )

df_long <- df %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("x"),
    names_to = "variable",
    values_to = "value"
  )

p_selection_overlay <- ggplot(df_long, aes(x = value, y = z)) +
  # draw "neither" first so colored points sit on top
  geom_point(
    data = dplyr::filter(df_long, grp == "neither"),
    color = "grey60", alpha = 0.6
  ) +
  geom_point(
    data = dplyr::filter(df_long, grp != "neither"),
    aes(color = grp), alpha = 0.9
  ) +
  scale_color_manual(values = c(
    both   = "green",
    BH_only = "blue",
    chai_only = "red"
  )) +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~variable, scales = "free_x") +
  theme_minimal() +
  labs(title = paste("Relationship between z and predictors", run_label))
save_plot("16s_wilcoxz_LLM_selection_overlay.pdf", p_selection_overlay, width = 10, height = 5)

p_space <- ggplot(df, aes(x = xA, y = xB, color = grp)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c(
      chai_only = "red",
      BH_only = "blue",
      both = "green",
      neither = "black"
    ),
    breaks = c("chai_only", "BH_only", "both", "neither"),
    labels = c("Only in chai", "Only in BH", "In both", "Neither")
  ) +
  theme_minimal() +
  labs(
    title = paste("LLM Space Map", run_label),
    subtitle = "Points colored by findings",
    x = "LLM 1",
    y = "LLM 2",
    color = "Group"
  ) +
  theme(legend.position = "right")
save_plot("16s_wilcoxz_LLM_space_map.pdf", p_space, width = 8, height = 6)

# Explore the selections
library(ape); library(phyloseq)

# inputs
tr <- phy_tree(genus)                 # phylo
sel_otus <- df_model$otu_rep[chai_sel]   # chai selected OTU IDs
# keep only OTUs that exist in the tree (safe)
sel_otus <- unique(intersect(sel_otus, tr$tip.label))

bh_otus <- df_model$otu_rep[bh_sel]   # BH selected OTU IDs

if (length(sel_otus) < 2) {
  warning("Too few selected OTUs mapped to tree (n = ", length(sel_otus),
          "). Skipping phylogenetic covariance heatmap.")
  R_phylo <- NULL
} else {
  # prune tree to selected OTUs
  tr_sel <- keep.tip(tr, sel_otus)

  # Ensure branch lengths exist for vcv.phylo()
  if (is.null(tr_sel$edge.length) || !length(tr_sel$edge.length) ||
      any(!is.finite(tr_sel$edge.length))) {
    warning("Selected tree has missing branch lengths; assigning Grafen branch lengths.")
    tr_sel <- ape::compute.brlen(tr_sel, method = "Grafen")
  }

  # Brownian-motion phylogenetic covariance, then correlation
  V <- vcv.phylo(tr_sel)     # covariance implied by tree branch lengths
  R_phylo <- cov2cor(V)      # correlation matrix

  # Rename OTU to taxa
  tax <- as.data.frame(tax_table(genus))
  tax$OTU <- rownames(tax)
  lab_df <- tax %>%
    filter(OTU %in% sel_otus) %>%
    mutate(
      raw = .data[["genus"]],
      clean = clean_taxon(raw),
      label = ifelse(is.na(clean) | clean == "", OTU, paste0(clean, " (", OTU, ")"))
    )

  lbl <- setNames(make.unique(lab_df$label), lab_df$OTU)
  rn <- rownames(R_phylo)
  cn <- colnames(R_phylo)
  rownames(R_phylo) <- ifelse(!is.na(lbl[rn]), lbl[rn], rn)
  colnames(R_phylo) <- ifelse(!is.na(lbl[cn]), lbl[cn], cn)

  library(pheatmap)
  pdf(plot_file("16s_wilcoxz_LLM_phylo_correlation_heatmap.pdf"), width = 10, height = 10)
  pheatmap(
    R_phylo,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    main = paste("Phylogenetic correlation among selected taxa", run_label)
  )
  dev.off()
}


# Tree
tr_all <- phy_tree(genus)

# map labels for ALL tips (not only selected)
tax_all <- as.data.frame(tax_table(genus))
tax_all$OTU <- rownames(tax_all)

rank_to_show <- "genus"

lab_all <- tax_all %>%
  mutate(
    raw = .data[[rank_to_show]],
    clean = clean_taxon(raw),
    label = ifelse(is.na(clean) | clean == "", OTU, clean)
  )

lbl_all <- setNames(make.unique(lab_all$label), lab_all$OTU)

# colors
tip_cols <- ifelse(tr_all$tip.label %in% sel_otus, "red", "black")

# plot with original OTU labels
pdf(plot_file("16s_wilcoxz_LLM_phylo_tree.pdf"), width = 12, height = 10)
plot(tr_all, show.tip.label = FALSE, cex = 0.7, tip.color = tip_cols,
     main = paste("Phylogenetic tree with selected genera by chai using LLM", run_label))

# add cleaned labels on top (so you see taxa names)
tiplabels(text = unname(lbl_all[tr_all$tip.label]) , tip = 1:length(tr_all$tip.label),
          frame = "none", adj = 0, cex = 0.7, col = tip_cols)
dev.off()


# Circular tree
library(dplyr)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

rank_to_show <- "genus"

anno <- tax_all %>%
  # tibble::rownames_to_column("OTU") %>%
  mutate(
    raw = .data[[rank_to_show]],
    clean = clean_taxon(raw),
    tip_label = ifelse(is.na(clean) | clean == "", OTU, clean),

    phylum_clean = clean_taxon(phylum),

    tip_color = dplyr::case_when(
      OTU %in% sel_otus & OTU %in% bh_otus ~ "both",
      OTU %in% sel_otus & !(OTU %in% bh_otus) ~ "chai_only",
      !(OTU %in% sel_otus) & OTU %in% bh_otus ~ "BH_only",
      TRUE ~ "neither"
    )  ) %>%
  dplyr::select(OTU, tip_label, phylum_clean, tip_color) %>%
  distinct(OTU, .keep_all = TRUE)

anno$tip_color <- factor(
  anno$tip_color,
  levels = c("chai_only", "BH_only", "both", "neither")
)

tr <- ape::as.phylo(tr_all)
p <- ggtree(tr, layout = "circular", size = 0.5)

p1 <- p %<+% anno

p2 <- p1 +
  new_scale_fill() +
  ggtreeExtra::geom_fruit(
    # REMOVE data = anno (it's already inside p2)
    geom = geom_tile,
    mapping = aes(fill = phylum_clean), # y = OTU is often inferred, but you can keep it if needed
    width = 0.2,
    offset = 0.05
  ) +
  labs(fill = "phylum")

p3 <- p2 +
  geom_tiplab(
    aes(label = tip_label, color = tip_color),
    size = 2,
    offset = 0.02,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name   = "Tip label group",
    values = c(chai_only="red", BH_only="blue", both="green", neither="black"),
    breaks = c("chai_only", "BH_only", "both", "neither"),
    labels = c("Only in chai", "Only in BH", "In both", "Neither"),
    drop   = FALSE
  ) +
  labs(title = paste("Phylogenetic tree with selected genera by chai on wilcoxon z with LLM", run_label)) +
  theme(legend.position = "right")

legend_df <- data.frame(
  x = NA_real_, y = NA_real_,
  tip_color = factor(c("chai_only","BH_only","both","neither"),
                     levels = c("chai_only","BH_only","both","neither"))
)

p4 <- p3 +
  geom_point(
    data = legend_df,
    aes(x = x, y = y, color = tip_color),
    inherit.aes = FALSE,
    size = 3
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = "right")
save_plot("16s_wilcoxz_LLM_phylo_tree_circular.pdf", p4, width = 12, height = 12)


################# Check which features have higher influence ###########################
full <- lab_all %>%
  left_join(gen_ag_gen, by = c("label" = "genus"))

top_sig_A <- full %>%
  dplyr::arrange(desc(abs(xA))) %>%
  head(20)
print(top_sig_A)
readr::write_csv(top_sig_A, plot_file("top_signal_A.csv"))

top_sig_B <- full %>%
  dplyr::arrange(desc(abs(xB))) %>%
  head(20)
print(top_sig_B)
readr::write_csv(top_sig_B, plot_file("top_signal_B.csv"))

top_sig_C <- full %>%
  dplyr::arrange(desc(abs(xC))) %>%
  head(20)
print(top_sig_C)
readr::write_csv(top_sig_C, plot_file("top_signal_C.csv"))
