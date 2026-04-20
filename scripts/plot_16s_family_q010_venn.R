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
  library(DESeq2)
  library(dplyr)
  library(VennDiagram)
  library(grid)
  library(locfdr)
  library(admix)
  library(mvtnorm)
  library(mixtools)
  library(KScorrect)
  library(mclust)
})

source(file.path(project_root, "Gopalakrishnan", "code", "conditionalParam.R"))
source(file.path(project_root, "Gopalakrishnan", "code", "naiveRemoveOneObs.R"))
source(file.path(project_root, "Gopalakrishnan", "code", "rGaussianMix.R"))
source(file.path(project_root, "Gopalakrishnan", "code", "utils.R"))
source(file.path(project_root, "Gopalakrishnan", "code", "chai.R"))
source(file.path(project_root, "Gopalakrishnan", "code", "clFDRselect.R"))

clean_taxon <- function(x) {
  y <- as.character(x)
  y <- sub("^.*?__", "", y)
  y <- gsub("[", "", y, fixed = TRUE)
  y <- gsub("]", "", y, fixed = TRUE)
  trimws(y)
}

output_dir <- file.path(project_root, "plots", "Gopalakrishnan")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

build_family_phyloseq <- function(project_root) {
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
  taxonomy <- taxonomy[otus_keep, , drop = FALSE]
  taxonomy <- taxonomy[match(tree$tip.label, rownames(taxonomy)), , drop = FALSE]
  otutable <- otutable[otus_keep, , drop = FALSE]

  ps <- phyloseq(
    otu_table(as.matrix(otutable), taxa_are_rows = TRUE),
    tax_table(as.matrix(taxonomy)),
    sample_data(meta_data),
    phy_tree(tree)
  )

  tax_glom(ps, taxrank = "family")
}

q_tag <- function(q) {
  sprintf("q%03d", as.integer(round(q * 100)))
}

compute_deseq2_pcoa <- function(family_ps, q_values = c(0.05, 0.10)) {
  family_pcoa <- cmdscale(cophenetic(phy_tree(family_ps)), k = 2)

  cts <- as(otu_table(family_ps), "matrix")
  md <- as(sample_data(family_ps), "data.frame")
  md$BinOutcomes <- factor(md$BinOutcomes, levels = c("NR", "R"))

  set.seed(123)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = md,
    design = ~ BinOutcomes
  )
  dds_res <- DESeq(dds, quiet = TRUE)
  res_05 <- results(dds_res, alpha = 0.05)

  set.seed(123)
  chai_fit <- chai(res_05$stat, family_pcoa, B = 1000)
  lfdr <- if (!is.null(chai_fit$lFDR)) chai_fit$lFDR else chai_fit$clFDR

  tax_df <- as.data.frame(tax_table(family_ps), stringsAsFactors = FALSE)
  tax_df$feature_id <- rownames(tax_df)
  tax_df$family_clean <- clean_taxon(tax_df$family)
  tax_df$clFDR <- as.numeric(lfdr)

  out <- vector("list", length(q_values))
  names(out) <- q_tag(q_values)

  for (i in seq_along(q_values)) {
    q <- q_values[[i]]
    sel <- clfdrselect(lfdr, q = q, max_clFDR = 1)
    selected <- tax_df[sel, c("feature_id", "kingdom", "phylum", "class", "order", "family", "family_clean", "clFDR")]
    selected <- selected[order(selected$clFDR), , drop = FALSE]

    out[[i]] <- list(
      selected_table = selected,
      selected_families = unique(selected$family_clean)
    )
  }

  out
}

load_wilcox_pcoa <- function(project_root, q) {
  path <- file.path(
    project_root,
    "plots",
    "Gopalakrishnan",
    sprintf("chai_selected_family_pcoa_%s.csv", q_tag(q))
  )
  if (!file.exists(path)) {
    stop("Expected Wilcoxon + PCoA selection file not found: ", path)
  }

  df <- read.csv(path, stringsAsFactors = FALSE)
  if (!"family_clean" %in% names(df)) {
    stop("Missing family_clean column in: ", path)
  }

  unique(df$family_clean)
}

load_legacy_llm <- function(project_root, q) {
  path <- file.path(
    project_root,
    "Gopalakrishnan",
    "plots",
    "featureSelectionLines",
    "combined",
    "weighted",
    sprintf("chai_selected_family_%s_aux_family_weighted_joined.csv", q_tag(q))
  )
  if (!file.exists(path)) {
    stop("Expected legacy LLM selection file not found: ", path)
  }

  df <- read.csv(path, stringsAsFactors = FALSE)
  if (!"family" %in% names(df)) {
    stop("Missing family column in: ", path)
  }

  unique(df$family)
}

family_ps <- build_family_phyloseq(project_root)
q_values <- c(0.10, 0.05)
deseq2_pcoa_by_q <- compute_deseq2_pcoa(family_ps, q_values = q_values)

for (q in q_values) {
  tag <- q_tag(q)
  deseq2_pcoa <- deseq2_pcoa_by_q[[tag]]
  wilcox_pcoa_families <- load_wilcox_pcoa(project_root, q)
  llm_wilcox_legacy_families <- load_legacy_llm(project_root, q)

  venn_sets <- list(
    "PCoA +\nWilcoxon" = sort(unique(wilcox_pcoa_families)),
    "PCoA +\nDESeq2" = sort(unique(deseq2_pcoa$selected_families)),
    "LLM +\nWilcoxon" = sort(unique(llm_wilcox_legacy_families))
  )

  membership_df <- data.frame(
    family = sort(unique(unlist(venn_sets))),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  for (set_name in names(venn_sets)) {
    membership_df[[set_name]] <- membership_df$family %in% venn_sets[[set_name]]
  }

  deseq2_out_csv <- file.path(output_dir, sprintf("chai_selected_family_deseq2_pcoa_%s.csv", tag))
  membership_out_csv <- file.path(output_dir, sprintf("16s_family_%s_venn_membership.csv", tag))
  plot_out_pdf <- file.path(output_dir, sprintf("16s_family_%s_venn.pdf", tag))

  write.csv(deseq2_pcoa$selected_table, deseq2_out_csv, row.names = FALSE)
  write.csv(membership_df, membership_out_csv, row.names = FALSE)

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
    main = sprintf("Selected families at q = %.2f", q),
    main.cex = 1.4,
    main.fontface = "bold"
  )

  pdf(plot_out_pdf, width = 8.5, height = 7)
  grid.newpage()
  grid.draw(venn_grob)
  dev.off()

  cat("Wrote Venn PDF:", plot_out_pdf, "\n")
  cat("Wrote DESeq2 + PCoA selected families CSV:", deseq2_out_csv, "\n")
  cat("Wrote membership CSV:", membership_out_csv, "\n")
  cat("Set sizes for", tag, ":\n")
  for (set_name in names(venn_sets)) {
    cat(" -", set_name, ":", length(venn_sets[[set_name]]), "\n")
  }
}
