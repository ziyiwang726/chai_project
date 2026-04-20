suppressPackageStartupMessages({
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

  library(ape)
  library(CATMicrobiome)
  library(phyloseq)
  library(DESeq2)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
})

output_dir <- file.path(project_root, "plots", "Gopalakrishnan")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    str_replace("^.*?__", "") %>%
    str_replace_all("\\[|\\]", "") %>%
    str_squish()
}

tip_group_levels <- c("chai_only", "BH_only", "both", "neither")
tip_group_labels <- c(
  chai_only = "Only in chai",
  BH_only = "Only in BH",
  both = "In both",
  neither = "Neither"
)
tip_group_palette <- c(
  chai_only = "#E69F00",
  BH_only = "#0072B2",
  both = "#CC79A7",
  neither = "#666666"
)
base_tip_label_cex <- 1.15
circular_tip_label_size <- 3.6
base_legend_cex <- 1.05
circular_legend_text_size <- 13
circular_legend_title_size <- 14

saved_chai_csv <- file.path(output_dir, "chai_selected_family_deseq2_pcoa_q005.csv")
saved_summary_rdata <- file.path(output_dir, "16s_DES_pcoa_family_updated_plot_data.RData")

if (!file.exists(saved_chai_csv)) {
  stop("Saved chai selection not found: ", saved_chai_csv)
}
if (!file.exists(saved_summary_rdata)) {
  stop("Saved summary RData not found: ", saved_summary_rdata)
}

otu_path <- system.file("extdata", "d1OTUtable.csv", package = "CATMicrobiome")
meta_path <- system.file("extdata", "d1Meta.csv", package = "CATMicrobiome")
tree_path <- system.file("extdata", "d1Tree.tree", package = "CATMicrobiome")
taxonomy_path <- file.path(project_root, "Gopalakrishnan", "LLMCode", "d1Taxonomy.csv")

if (!file.exists(taxonomy_path)) {
  stop("Expected taxonomy CSV not found: ", taxonomy_path)
}

otutable <- read.csv(otu_path, header = TRUE, row.names = 1)
meta_data <- read.csv(meta_path, header = TRUE, row.names = 1)
taxonomy <- read.csv(taxonomy_path, header = TRUE, row.names = 1)
tree <- read.tree(tree_path)

otus_keep <- intersect(tree$tip.label, rownames(taxonomy))
tree <- keep.tip(tree, otus_keep)
taxonomy <- taxonomy[match(tree$tip.label, rownames(taxonomy)), , drop = FALSE]
otutable <- otutable[otus_keep, , drop = FALSE]

ps <- phyloseq(
  otu_table(as.matrix(otutable), taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy)),
  sample_data(meta_data),
  phy_tree(tree)
)
family <- tax_glom(ps, taxrank = "family")

cts <- as(otu_table(family), "matrix")
md <- as(sample_data(family), "data.frame")
md$BinOutcomes <- factor(md$BinOutcomes, levels = c("NR", "R"))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = md, design = ~ BinOutcomes)
dds_res <- DESeq(dds, quiet = TRUE)
res_05 <- results(dds_res, alpha = 0.05)
des_p_adjusted <- p.adjust(res_05$pvalue, method = "BH")
sel_otus_bh <- rownames(res_05)[which(des_p_adjusted <= 0.05)]

chai_selected <- read.csv(saved_chai_csv, stringsAsFactors = FALSE)
sel_otus_chai <- chai_selected$feature_id

load(saved_summary_rdata)
q005_summary <- gopa_sum_16s_DES[gopa_sum_16s_DES$q == 0.05, c("chai", "BH")]
stopifnot(length(sel_otus_chai) == q005_summary$chai)
stopifnot(length(sel_otus_bh) == q005_summary$BH)

tr_all <- phy_tree(family)
sel_otus_chai <- intersect(sel_otus_chai, tr_all$tip.label)
sel_otus_bh <- intersect(sel_otus_bh, tr_all$tip.label)

tax_all <- as.data.frame(tax_table(family), stringsAsFactors = FALSE)
tax_all$OTU <- rownames(tax_all)

lab_all <- tax_all %>%
  mutate(
    clean = clean_taxon(family),
    label = ifelse(is.na(clean) | clean == "", OTU, clean)
  )
lbl_all <- setNames(make.unique(lab_all$label), lab_all$OTU)

tip_group <- dplyr::case_when(
  tr_all$tip.label %in% sel_otus_chai & tr_all$tip.label %in% sel_otus_bh ~ "both",
  tr_all$tip.label %in% sel_otus_chai & !(tr_all$tip.label %in% sel_otus_bh) ~ "chai_only",
  !(tr_all$tip.label %in% sel_otus_chai) & tr_all$tip.label %in% sel_otus_bh ~ "BH_only",
  TRUE ~ "neither"
)
tip_cols <- unname(tip_group_palette[tip_group])

tip_group_df <- tax_all %>%
  transmute(
    feature_id = OTU,
    family_clean = clean_taxon(family),
    phylum_clean = clean_taxon(phylum),
    tip_group = factor(
      dplyr::case_when(
        OTU %in% sel_otus_chai & OTU %in% sel_otus_bh ~ "both",
        OTU %in% sel_otus_chai & !(OTU %in% sel_otus_bh) ~ "chai_only",
        !(OTU %in% sel_otus_chai) & OTU %in% sel_otus_bh ~ "BH_only",
        TRUE ~ "neither"
      ),
      levels = tip_group_levels
    )
  )
write.csv(
  tip_group_df,
  file.path(output_dir, "16s_DES_BH_PCoA_phylo_tree_tip_groups.csv"),
  row.names = FALSE
)

pdf(file.path(output_dir, "16s_DES_BH_PCoA_phylo_tree.pdf"), width = 12, height = 10)
plot(
  tr_all,
  show.tip.label = FALSE,
  cex = 0.7,
  tip.color = tip_cols,
  main = "Phylogenetic tree with selected families by chai on DESeq2 z with PCoA (BH, q = 0.05)"
)
tiplabels(
  text = unname(lbl_all[tr_all$tip.label]),
  tip = seq_along(tr_all$tip.label),
  frame = "none",
  adj = 0,
  cex = base_tip_label_cex,
  col = tip_cols
)
legend(
  "bottomright",
  legend = unname(tip_group_labels[tip_group_levels]),
  col = unname(tip_group_palette[tip_group_levels]),
  pch = 16,
  bty = "n",
  cex = base_legend_cex,
  pt.cex = 0.9,
  x.intersp = 0.5,
  y.intersp = 0.2,
  inset = 0.02
)
dev.off()

anno_des <- tax_all %>%
  mutate(
    tip_label = ifelse(is.na(clean_taxon(family)) | clean_taxon(family) == "", OTU, clean_taxon(family)),
    phylum_clean = clean_taxon(phylum),
    tip_color = factor(
      dplyr::case_when(
        OTU %in% sel_otus_chai & OTU %in% sel_otus_bh ~ "both",
        OTU %in% sel_otus_chai & !(OTU %in% sel_otus_bh) ~ "chai_only",
        !(OTU %in% sel_otus_chai) & OTU %in% sel_otus_bh ~ "BH_only",
        TRUE ~ "neither"
      ),
      levels = tip_group_levels
    )
  ) %>%
  dplyr::select(OTU, tip_label, phylum_clean, tip_color) %>%
  distinct(OTU, .keep_all = TRUE)

tr <- ape::as.phylo(tr_all)
p <- ggtree(tr, layout = "circular", size = 0.5)

p_circular <- p %<+% anno_des +
  new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom = geom_tile,
    mapping = aes(fill = phylum_clean),
    width = 0.2,
    offset = 0.05
  ) +
  labs(fill = "phylum") +
  geom_tiplab(
    aes(label = tip_label, color = tip_color),
    size = circular_tip_label_size,
    offset = 0.02,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name = "Tip label group",
    values = tip_group_palette,
    breaks = tip_group_levels,
    labels = unname(tip_group_labels[tip_group_levels]),
    drop = FALSE
  ) +
  labs(title = "Phylogenetic tree with selected families by chai on DESeq2 z with PCoA (BH, q = 0.05)") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = circular_legend_text_size),
    legend.title = element_text(size = circular_legend_title_size)
  )

legend_df <- data.frame(
  x = NA_real_,
  y = NA_real_,
  tip_color = factor(tip_group_levels, levels = tip_group_levels)
)

p_circular <- p_circular +
  geom_point(
    data = legend_df,
    aes(x = x, y = y, color = tip_color),
    inherit.aes = FALSE,
    size = 3
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = circular_legend_text_size),
    legend.title = element_text(size = circular_legend_title_size)
  )

ggsave(
  file.path(output_dir, "16s_DES_BH_PCoA_phylo_tree_circular.pdf"),
  p_circular,
  width = 12,
  height = 12
)

message("Saved tree plots to ", output_dir)
