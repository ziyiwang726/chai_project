#!/usr/bin/env Rscript

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  idx <- grep("^--file=", cmd_args)
  if (length(idx)) return(normalizePath(sub("^--file=", "", cmd_args[idx[1]])))
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(normalizePath(ofile))
  NA_character_
}

script_path <- get_script_path()
project_root <- if (!is.na(script_path)) dirname(dirname(script_path)) else normalizePath(getwd())
project_r_libs <- file.path(project_root, "r_libs")
if (dir.exists(project_r_libs)) {
  .libPaths(c(normalizePath(project_r_libs), .libPaths()))
}

suppressPackageStartupMessages({
  library(ape)
  library(CATMicrobiome)
  library(phyloseq)
  library(coin)
  library(DESeq2)
  library(locfdr)
  library(admix)
  library(mvtnorm)
  library(mixtools)
  library(KScorrect)
  library(mclust)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(VennDiagram)
  library(grid)
  library(adaptMT)
  library(AdaPTGMM)
  library(IHW)
  library(splines)
})

source(file.path(project_root, "Gopalakrishnan", "code", "conditionalParam.R"))
source(file.path(project_root, "Gopalakrishnan", "code", "naiveRemoveOneObs.R"))
source(file.path(project_root, "Gopalakrishnan", "code", "rGaussianMix.R"))
source(file.path(project_root, "Gopalakrishnan", "code", "utils.R"))
source(file.path(project_root, "Gopalakrishnan", "code", "chai.R"))
source(file.path(project_root, "Gopalakrishnan", "code", "clFDRselect.R"))
source(file.path(project_root, "Gopalakrishnan", "code", "color_helper.R"))

if (!exists("lFDRselect")) {
  lFDRselect <- clfdrselect
}

output_dir <- file.path(project_root, "plots", "Gopalakrishnan")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    str_replace("^.*?__", "") %>%
    str_replace_all("\\[|\\]", "") %>%
    str_squish()
}

q_tag <- function(q) {
  sprintf("q%03d", as.integer(round(q * 100)))
}

build_ps <- function(project_root) {
  otu_path <- system.file("extdata", "d1OTUtable.csv", package = "CATMicrobiome")
  meta_path <- system.file("extdata", "d1Meta.csv", package = "CATMicrobiome")
  tree_path <- system.file("extdata", "d1Tree.tree", package = "CATMicrobiome")
  taxonomy_path <- file.path(project_root, "Gopalakrishnan", "LLMCode", "d1Taxonomy.csv")

  if (!file.exists(taxonomy_path)) {
    stop("Expected taxonomy CSV not found: ", taxonomy_path)
  }

  otutable <- read.csv(otu_path, header = TRUE, row.names = 1)
  taxonomy <- read.csv(taxonomy_path, header = TRUE, row.names = 1)
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

compute_wilcox_pz <- function(t_abun, outcome) {
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

make_rank_table <- function(rank_ps, coords, p_values, z_values, lfdr, rank_col) {
  ord <- order(lfdr)
  avg_fdr <- cumsum(lfdr[ord]) / seq_along(ord)
  avg_lookup <- rep(NA_real_, length(lfdr))
  avg_lookup[ord] <- avg_fdr
  rank_lookup <- integer(length(lfdr))
  rank_lookup[ord] <- seq_along(ord)

  tax_df <- as.data.frame(tax_table(rank_ps), stringsAsFactors = FALSE)
  tax_df$feature_id <- rownames(tax_df)
  tax_df$rank_clean <- clean_taxon(tax_df[[rank_col]])
  tax_df$PCoA1 <- coords[, 1]
  tax_df$PCoA2 <- coords[, 2]
  tax_df$p_value <- as.numeric(p_values)
  tax_df$z_value <- as.numeric(z_values)
  tax_df$clFDR <- as.numeric(lfdr)
  tax_df$clFDR_rank <- rank_lookup
  tax_df$avgFDR_at_rank <- avg_lookup
  tax_df %>% arrange(clFDR)
}

stouffer_nonzero <- function(v) {
  v2 <- v[!is.na(v) & v != 0]
  if (!length(v2)) return(0)
  stouffer_z <- function(z, ignore.na = TRUE) {
    if (ignore.na) {
      z <- z[!is.na(z)]
    }
    if (!length(z)) {
      return(0)
    }
    sum(z) / sqrt(length(z))
  }
  stouffer_z(v2, ignore.na = FALSE)
}

make_aux_tax_level <- function(aux_otu, level, aux_covariates) {
  aux_otu %>%
    filter(!is.na(.data[[level]]), .data[[level]] != "") %>%
    group_by(.data[[level]]) %>%
    summarise(
      across(all_of(aux_covariates), stouffer_nonzero, .names = "x{col}"),
      .groups = "drop"
    ) %>%
    mutate(taxon_key = clean_taxon(.data[[level]]))
}

compute_ihw_rejections <- function(p, covariate, q_levels) {
  finite_vals <- covariate[is.finite(covariate)]
  if (length(unique(finite_vals)) <= 1) {
    return(rep(NA_real_, length(q_levels)))
  }

  sapply(q_levels, function(q) {
    ihw_fit <- ihw(pvalues = p, covariates = covariate, alpha = q, nbins = 5)
    rejections(ihw_fit)
  })
}

fit_ordershape <- function(p, covariate) {
  if (!requireNamespace("OrderShapeEM", quietly = TRUE)) {
    return(NULL)
  }
  finite_vals <- covariate[is.finite(covariate)]
  if (length(unique(finite_vals)) <= 1) {
    return(NULL)
  }

  OrderShapeEM::OrderShapeEM(
    pvals = p,
    order.var = covariate,
    OrderShapeEM::OrderShapeEM.control(trace = TRUE)
  )
}

extract_selected_taxa <- function(selection_table, column_name) {
  sort(unique(selection_table[[column_name]]))
}

ps <- build_ps(project_root)
rank_name <- "genus"
rank_ps <- tax_glom(ps, taxrank = rank_name, NArm = FALSE)
rank_dist <- cophenetic(phy_tree(rank_ps))
rank_pcoa <- cmdscale(rank_dist, k = 2)
colnames(rank_pcoa) <- c("X1", "X2")

cs <- colSums(otu_table(rank_ps))
cs[cs == 0] <- 1
rank_abun <- sweep(otu_table(rank_ps), 2, cs, "/") * 100
t_abun <- t(rank_abun)
outcome <- factor(sample_data(rank_ps)$BinOutcomes, levels = c("NR", "R"))

wilcox_pz <- compute_wilcox_pz(t_abun, outcome)
wilcox_p <- wilcox_pz$p_value
wilcox_z <- wilcox_pz$z_value

set.seed(123)
chai_wilcox_pcoa <- chai(wilcox_z, rank_pcoa, B = 1000)
lfdr_wilcox_pcoa <- if (!is.null(chai_wilcox_pcoa$lFDR)) chai_wilcox_pcoa$lFDR else chai_wilcox_pcoa$clFDR

genus_pcoa_table <- make_rank_table(
  rank_ps = rank_ps,
  coords = rank_pcoa,
  p_values = wilcox_p,
  z_values = wilcox_z,
  lfdr = lfdr_wilcox_pcoa,
  rank_col = rank_name
)

wilcox_sel_q05 <- lFDRselect(lfdr_wilcox_pcoa, 0.05, 1)
wilcox_sel_q10 <- lFDRselect(lfdr_wilcox_pcoa, 0.10, 1)
chai_selected_genus_pcoa_q005 <- genus_pcoa_table[wilcox_sel_q05, , drop = FALSE] %>% arrange(clFDR)
chai_selected_genus_pcoa_q010 <- genus_pcoa_table[wilcox_sel_q10, , drop = FALSE] %>% arrange(clFDR)

write.csv(
  chai_selected_genus_pcoa_q005,
  file.path(output_dir, "chai_selected_genus_pcoa_q005.csv"),
  row.names = FALSE
)
write.csv(
  chai_selected_genus_pcoa_q010,
  file.path(output_dir, "chai_selected_genus_pcoa_q010.csv"),
  row.names = FALSE
)

q_levels <- seq(0.01, 0.10, by = 0.01)
x <- data.frame(X1 = rank_pcoa[, 1], X2 = rank_pcoa[, 2])
formula_ns <- list(
  "~ns(X1, df = 2) + ns(X2, df = 2)",
  "~ns(X1, df = 4) + ns(X2, df = 4)",
  "~ns(X1, df = 6) + ns(X2, df = 6)"
)

wilcox_p_adjusted <- p.adjust(wilcox_p, method = "BH")

glm_wilcox_ns <- adapt_glm(x, pvals = wilcox_p, alphas = q_levels, pi_formulas = formula_ns, mu_formulas = formula_ns)
gmm_wilcox_ns_p <- adapt_gmm(x, pvals = wilcox_p, alphas = q_levels, beta_formulas = formula_ns)
gmm_wilcox_ns_z <- adapt_gmm(x, z = wilcox_z, alphas = q_levels, beta_formulas = formula_ns, testing = "two_sided")

ordershape_wilcox_x1 <- fit_ordershape(wilcox_p, rank_pcoa[, 1])
ordershape_wilcox_x2 <- fit_ordershape(wilcox_p, rank_pcoa[, 2])

fdr_wilcox_theo <- if (requireNamespace("FDRreg", quietly = TRUE)) {
  set.seed(123)
  FDRreg::FDRreg(wilcox_z, rank_pcoa, nulltype = "theoretical")
} else {
  NULL
}

chai_wilcox_rejs <- sapply(q_levels, function(q) length(lFDRselect(lfdr_wilcox_pcoa, q, 1)))
adapt_glm_wilcox_rej <- glm_wilcox_ns$nrejs[seq_along(q_levels)]
adapt_gmm_wilcox_p_rej <- gmm_wilcox_ns_p$nrejs[seq_along(q_levels)]
adapt_gmm_wilcox_z_rej <- gmm_wilcox_ns_z$nrejs[seq_along(q_levels)]
ordershape_wilcox_x1_rej <- if (!is.null(ordershape_wilcox_x1)) sapply(q_levels, function(q) sum(ordershape_wilcox_x1$fdr <= q)) else rep(NA_real_, length(q_levels))
ordershape_wilcox_x2_rej <- if (!is.null(ordershape_wilcox_x2)) sapply(q_levels, function(q) sum(ordershape_wilcox_x2$fdr <= q)) else rep(NA_real_, length(q_levels))
fdrreg_wilcox_rej <- if (!is.null(fdr_wilcox_theo)) sapply(q_levels, function(q) length(which(fdr_wilcox_theo$FDR <= q))) else rep(NA_real_, length(q_levels))
ihw_wilcox_x1_rej <- compute_ihw_rejections(wilcox_p, rank_pcoa[, 1], q_levels)
ihw_wilcox_x2_rej <- compute_ihw_rejections(wilcox_p, rank_pcoa[, 2], q_levels)
bh_wilcox_rej <- sapply(q_levels, function(q) length(which(wilcox_p_adjusted <= q)))

gopa_sum_16s_genus <- data.frame(
  q = q_levels,
  chai = chai_wilcox_rejs,
  adapt_glm = adapt_glm_wilcox_rej,
  adapt_gmm_p = adapt_gmm_wilcox_p_rej,
  adapt_gmm_z = adapt_gmm_wilcox_z_rej,
  OrderShapeEM_x1 = ordershape_wilcox_x1_rej,
  OrderShapeEM_x2 = ordershape_wilcox_x2_rej,
  FDRreg = fdrreg_wilcox_rej,
  IHW_x1 = ihw_wilcox_x1_rej,
  IHW_x2 = ihw_wilcox_x2_rej,
  BH = bh_wilcox_rej
)

gopa_16s_genus_long <- gopa_sum_16s_genus %>%
  pivot_longer(cols = -q, names_to = "Method", values_to = "Rejections") %>%
  mutate(Method = as.character(Method))

p_gopa_16s_pcoa_genus <- plot_rejections_vs_q(
  gopa_16s_genus_long,
  title = "Number of Discoveries vs q, on Wilcoxon z-statistics",
  x_breaks = sort(unique(gopa_16s_genus_long$q))
)

pdf(file.path(output_dir, "16s_wilcoxz_pcoa_genus.pdf"), width = 8, height = 4)
plot(p_gopa_16s_pcoa_genus)
dev.off()

save(
  gopa_sum_16s_genus,
  gopa_16s_genus_long,
  p_gopa_16s_pcoa_genus,
  file = file.path(output_dir, "16s_wilcoxz_pcoa_genus_plot_data.RData"),
  compress = "xz"
)

cts <- as(otu_table(rank_ps), "matrix")
md <- as(sample_data(rank_ps), "data.frame")
md$BinOutcomes <- factor(md$BinOutcomes, levels = c("NR", "R"))

set.seed(123)
dds <- DESeqDataSetFromMatrix(countData = cts, colData = md, design = ~ BinOutcomes)
dds_res <- DESeq(dds, quiet = TRUE)
res_05 <- results(dds_res, alpha = 0.05)

des_stat <- res_05$stat
des_p <- res_05$pvalue
des_p_adjusted <- p.adjust(des_p, method = "BH")

set.seed(123)
chai_des_pcoa <- chai(des_stat, rank_pcoa, B = 1000)
lfdr_des_pcoa <- if (!is.null(chai_des_pcoa$lFDR)) chai_des_pcoa$lFDR else chai_des_pcoa$clFDR

genus_des_table <- as.data.frame(tax_table(rank_ps), stringsAsFactors = FALSE)
genus_des_table$feature_id <- rownames(genus_des_table)
genus_des_table$genus_clean <- clean_taxon(genus_des_table$genus)
genus_des_table$clFDR <- as.numeric(lfdr_des_pcoa)
genus_des_table <- genus_des_table %>% arrange(clFDR)

des_sel_q05 <- lFDRselect(lfdr_des_pcoa, 0.05, 1)
des_sel_q10 <- lFDRselect(lfdr_des_pcoa, 0.10, 1)
chai_selected_genus_deseq2_pcoa_q005 <- genus_des_table[des_sel_q05, c("feature_id", "kingdom", "phylum", "class", "order", "family", "genus", "genus_clean", "clFDR"), drop = FALSE] %>% arrange(clFDR)
chai_selected_genus_deseq2_pcoa_q010 <- genus_des_table[des_sel_q10, c("feature_id", "kingdom", "phylum", "class", "order", "family", "genus", "genus_clean", "clFDR"), drop = FALSE] %>% arrange(clFDR)

write.csv(
  chai_selected_genus_deseq2_pcoa_q005,
  file.path(output_dir, "chai_selected_genus_deseq2_pcoa_q005.csv"),
  row.names = FALSE
)
write.csv(
  chai_selected_genus_deseq2_pcoa_q010,
  file.path(output_dir, "chai_selected_genus_deseq2_pcoa_q010.csv"),
  row.names = FALSE
)

glm_des_ns <- adapt_glm(x, pvals = des_p, alphas = q_levels, pi_formulas = formula_ns, mu_formulas = formula_ns)
gmm_des_ns_p <- adapt_gmm(x, pvals = des_p, alphas = q_levels, beta_formulas = formula_ns)
gmm_des_ns_z <- adapt_gmm(x, z = des_stat, alphas = q_levels, beta_formulas = formula_ns, testing = "two_sided")

ordershape_des_x1 <- fit_ordershape(des_p, rank_pcoa[, 1])
ordershape_des_x2 <- fit_ordershape(des_p, rank_pcoa[, 2])

fdr_des_theo <- if (requireNamespace("FDRreg", quietly = TRUE)) {
  set.seed(123)
  FDRreg::FDRreg(des_stat, rank_pcoa, nulltype = "theoretical")
} else {
  NULL
}

chai_des_rejs <- sapply(q_levels, function(q) length(lFDRselect(lfdr_des_pcoa, q, 1)))
adapt_glm_des_rej <- glm_des_ns$nrejs[seq_along(q_levels)]
adapt_gmm_des_p_rej <- gmm_des_ns_p$nrejs[seq_along(q_levels)]
adapt_gmm_des_z_rej <- gmm_des_ns_z$nrejs[seq_along(q_levels)]
ordershape_des_x1_rej <- if (!is.null(ordershape_des_x1)) sapply(q_levels, function(q) sum(ordershape_des_x1$fdr <= q)) else rep(NA_real_, length(q_levels))
ordershape_des_x2_rej <- if (!is.null(ordershape_des_x2)) sapply(q_levels, function(q) sum(ordershape_des_x2$fdr <= q)) else rep(NA_real_, length(q_levels))
fdrreg_des_rej <- if (!is.null(fdr_des_theo)) sapply(q_levels, function(q) length(which(fdr_des_theo$FDR <= q))) else rep(NA_real_, length(q_levels))
ihw_des_x1_rej <- compute_ihw_rejections(des_p, rank_pcoa[, 1], q_levels)
ihw_des_x2_rej <- compute_ihw_rejections(des_p, rank_pcoa[, 2], q_levels)
bh_des_rej <- sapply(q_levels, function(q) length(which(des_p_adjusted <= q)))

gopa_sum_16s_DES_genus <- data.frame(
  q = q_levels,
  chai = chai_des_rejs,
  adapt_glm = adapt_glm_des_rej,
  adapt_gmm_p = adapt_gmm_des_p_rej,
  adapt_gmm_z = adapt_gmm_des_z_rej,
  OrderShapeEM_x1 = ordershape_des_x1_rej,
  OrderShapeEM_x2 = ordershape_des_x2_rej,
  FDRreg = fdrreg_des_rej,
  IHW_x1 = ihw_des_x1_rej,
  IHW_x2 = ihw_des_x2_rej,
  BH = bh_des_rej
)

gopa_16s_DES_genus_long <- gopa_sum_16s_DES_genus %>%
  pivot_longer(cols = -q, names_to = "Method", values_to = "Rejections") %>%
  mutate(Method = as.character(Method))

p_gopa_16s_pcoa_DES_genus <- plot_rejections_vs_q(
  gopa_16s_DES_genus_long,
  title = "Number of Discoveries vs q, on DESeq2 statistics",
  x_breaks = sort(unique(gopa_16s_DES_genus_long$q))
)

pdf(file.path(output_dir, "16s_DES_pcoa_genus_updated.pdf"), width = 8, height = 4)
plot(p_gopa_16s_pcoa_DES_genus)
dev.off()

save(
  gopa_sum_16s_DES_genus,
  gopa_16s_DES_genus_long,
  p_gopa_16s_pcoa_DES_genus,
  file = file.path(output_dir, "16s_DES_pcoa_genus_updated_plot_data.RData"),
  compress = "xz"
)

aux_file <- file.path(
  project_root,
  "Gopalakrishnan",
  "sideCov",
  "combined",
  "weighted",
  "selectionTargetGenus",
  "d1Taxonomy_genus_ABCD_genus_and_below_weighted.csv"
)
if (!file.exists(aux_file)) {
  stop("Expected weighted genus auxiliary CSV not found: ", aux_file)
}

main_df <- {
  tax_df <- as.data.frame(tax_table(rank_ps)) %>%
    tibble::rownames_to_column("otu_rep") %>%
    transmute(otu_rep, rank_value = as.character(genus))

  pz_map <- wilcox_pz %>%
    transmute(
      otu_rep = feature,
      p_main = as.numeric(p_value),
      z_main = as.numeric(z_value)
    )

  pz_map %>%
    left_join(tax_df, by = "otu_rep") %>%
    mutate(taxon_key = clean_taxon(rank_value))
}

aux_otu <- readr::read_csv(aux_file, show_col_types = FALSE)
aux_rank <- make_aux_tax_level(aux_otu, level = "genus", aux_covariates = c("A", "B", "C"))
df_model <- main_df %>%
  left_join(aux_rank %>% select(taxon_key, xA, xB, xC), by = "taxon_key") %>%
  mutate(across(c(xA, xB, xC), ~ tidyr::replace_na(.x, 0)))

X_aux <- as.data.frame(df_model[, c("xA", "xB", "xC"), drop = FALSE])
keep_aux <- complete.cases(df_model$z_main, X_aux)
df_model <- df_model[keep_aux, , drop = FALSE]
X_aux <- X_aux[keep_aux, , drop = FALSE]

set.seed(123)
chai_llm_weighted <- chai(df_model$z_main, X_aux, B = 1000)
lfdr_llm_weighted <- if (!is.null(chai_llm_weighted$lFDR)) chai_llm_weighted$lFDR else chai_llm_weighted$clFDR
llm_sel_q05 <- lFDRselect(lfdr_llm_weighted, 0.05, 1)
llm_sel_q10 <- lFDRselect(lfdr_llm_weighted, 0.10, 1)

llm_genus_q05 <- sort(unique(df_model$taxon_key[llm_sel_q05]))
llm_genus_q10 <- sort(unique(df_model$taxon_key[llm_sel_q10]))
wilcox_genus_q05 <- extract_selected_taxa(chai_selected_genus_pcoa_q005, "rank_clean")
wilcox_genus_q10 <- extract_selected_taxa(chai_selected_genus_pcoa_q010, "rank_clean")
deseq_genus_q05 <- extract_selected_taxa(chai_selected_genus_deseq2_pcoa_q005, "genus_clean")
deseq_genus_q10 <- extract_selected_taxa(chai_selected_genus_deseq2_pcoa_q010, "genus_clean")

venn_specs <- list(
  q005 = list(
    q = 0.05,
    pcoa_wilcox = wilcox_genus_q05,
    pcoa_deseq = deseq_genus_q05,
    llm_wilcox = llm_genus_q05
  ),
  q010 = list(
    q = 0.10,
    pcoa_wilcox = wilcox_genus_q10,
    pcoa_deseq = deseq_genus_q10,
    llm_wilcox = llm_genus_q10
  )
)

for (tag in names(venn_specs)) {
  spec <- venn_specs[[tag]]
  venn_sets <- list(
    "PCoA +\nWilcoxon" = spec$pcoa_wilcox,
    "PCoA +\nDESeq2" = spec$pcoa_deseq,
    "LLM +\nWilcoxon" = spec$llm_wilcox
  )

  membership_df <- data.frame(
    genus = sort(unique(unlist(venn_sets))),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  for (set_name in names(venn_sets)) {
    membership_df[[set_name]] <- membership_df$genus %in% venn_sets[[set_name]]
  }

  write.csv(
    membership_df,
    file.path(output_dir, sprintf("16s_genus_%s_venn_membership.csv", tag)),
    row.names = FALSE
  )

  venn_grob <- venn.diagram(
    x = venn_sets,
    category.names = names(venn_sets),
    filename = NULL,
    fill = c("#D55E00", "#009E73", "#56B4E9"),
    alpha = c(0.38, 0.38, 0.38),
    lwd = 2,
    col = c("#B24D00", "#007B5F", "#3A8CC1"),
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.1,
    cat.fontface = "bold",
    cat.fontfamily = "sans",
    cat.pos = c(-22, 22, 180),
    cat.dist = c(0.08, 0.08, 0.08),
    margin = 0.08,
    main = sprintf("Selected genera at q = %.2f", spec$q),
    main.cex = 1.4,
    main.fontface = "bold"
  )

  pdf(file.path(output_dir, sprintf("16s_genus_%s_venn.pdf", tag)), width = 8.5, height = 7)
  grid.newpage()
  grid.draw(venn_grob)
  dev.off()
}

cat("Wrote:", file.path(output_dir, "16s_DES_pcoa_genus_updated.pdf"), "\n")
cat("Wrote:", file.path(output_dir, "16s_genus_q005_venn.pdf"), "\n")
cat("Wrote:", file.path(output_dir, "16s_genus_q010_venn.pdf"), "\n")
