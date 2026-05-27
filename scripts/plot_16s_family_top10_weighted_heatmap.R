script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
project_root <- normalizePath(file.path(dirname(script_path), ".."))

project_lib <- file.path(project_root, "r_libs")
if (dir.exists(project_lib)) {
  .libPaths(c(normalizePath(project_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(ape)
  library(coin)
  library(dplyr)
  library(pheatmap)
  library(phyloseq)
  library(readr)
  library(stringr)
  library(tibble)
})

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

  p_z_values <- as.data.frame(t(p_z_values))
  p_z_values$feature <- rownames(p_z_values)
  rownames(p_z_values) <- NULL
  p_z_values
}

make_main_df <- function(t_abun, pz_df, ps_obj, rank) {
  tax_df <- as.data.frame(tax_table(ps_obj)) %>%
    rownames_to_column("otu_rep")

  rank_col <- names(tax_df)[tolower(names(tax_df)) == tolower(rank)]
  if (length(rank_col) != 1) {
    stop(
      "Could not uniquely find taxonomy column for rank = '", rank, "'. ",
      "Available columns: ", paste(names(tax_df), collapse = ", ")
    )
  }

  tax_map <- tax_df %>%
    transmute(
      otu_rep,
      rank_value = as.character(.data[[rank_col]])
    )

  pz_map <- pz_df %>%
    transmute(
      otu_rep = feature,
      p_main = as.numeric(p_value),
      z_main = as.numeric(z_value)
    )

  pz_map %>%
    left_join(tax_map, by = "otu_rep") %>%
    mutate(taxon_key = clean_taxon(rank_value))
}

otu_path <- system.file("extdata", "d1OTUtable.csv", package = "CATMicrobiome")
taxonomy_path <- system.file("extdata", "d1Taxonomy.csv", package = "CATMicrobiome")
meta_path <- system.file("extdata", "d1Meta.csv", package = "CATMicrobiome")
tree_path <- system.file("extdata", "d1Tree.tree", package = "CATMicrobiome")

otutable <- read.csv(otu_path, header = TRUE, row.names = 1)
taxonomy <- read.csv(taxonomy_path, header = TRUE, row.names = 1)
meta_data <- read.csv(meta_path, header = TRUE, row.names = 1)
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

family_abun <- otu_table(family)
cs <- colSums(family_abun)
cs[cs == 0] <- 1
family_abun <- sweep(family_abun, 2, cs, "/") * 100
t_family_abun <- t(family_abun)

outcome <- factor(sample_data(family)$BinOutcomes, levels = c("NR", "R"))
family_pz <- compute_pz_wilcox(t_family_abun, outcome)
family_pz_tax <- make_main_df(t_family_abun, family_pz, family, rank = "family")

side_info_path <- file.path(
  project_root,
  "Gopalakrishnan",
  "sideCov",
  "combined",
  "weighted",
  "selectionTargetFamily",
  "d1Taxonomy_family_ABCD_family_and_below_weighted.csv"
)

side_info <- read_csv(side_info_path, show_col_types = FALSE) %>%
  transmute(
    taxon_key = clean_taxon(family),
    A = as.numeric(A),
    B = as.numeric(B),
    C = as.numeric(C)
  ) %>%
  distinct(taxon_key, .keep_all = TRUE)

eligible_df <- family_pz_tax %>%
  left_join(side_info, by = "taxon_key")

ranked_df <- eligible_df %>%
  filter(
    !is.na(z_main),
    !is.na(A),
    !is.na(B),
    !is.na(C)
  ) %>%
  mutate(abs_z = abs(z_main)) %>%
  arrange(desc(abs_z), desc(z_main), taxon_key)

if (nrow(ranked_df) < 20) {
  stop("Only found ", nrow(ranked_df), " family-level features with non-missing z and LLM side information.")
}

global_mat <- log10(t_family_abun[, ranked_df$otu_rep, drop = FALSE] + 1e-6)
body_limits <- range(global_mat, finite = TRUE)
body_palette <- colorRampPalette(c("#FFFFFF", "#FEE8C8", "#FC8D59", "#B30000"))(100)
body_breaks <- seq(floor(body_limits[1]), ceiling(body_limits[2]), length.out = length(body_palette) + 1)
legend_breaks <- pretty(body_limits, n = 5)

output_dir <- file.path(project_root, "plots", "Gopalakrishnan")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

annotation_colors <- list(
  Group = c("NR" = "#33A02C", "R" = "#6A3D9A"),
  Direction = c("Higher in R" = "#D73027", "Higher in NR" = "#4575B4"),
  "Auxiliary information" = colorRampPalette(c("#2C7FB8", "#F7F7F7", "#FDB863"))(100)
)

write_padded_heatmap <- function(gt, pdf_out, png_out, width_in, height_in) {
  grDevices::pdf(pdf_out, width = width_in, height = height_in, useDingbats = FALSE)
  grid::grid.newpage()
  grid::grid.draw(gt)
  grDevices::dev.off()

  grDevices::png(png_out, width = width_in, height = height_in, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(gt)
  grDevices::dev.off()
}

write_heatmap_set <- function(selected_df, slug, title_text) {
  abun_sel <- t_family_abun[, selected_df$otu_rep, drop = FALSE]
  colnames(abun_sel) <- selected_df$taxon_key

  row_ord <- order(outcome)
  outcome_ord <- outcome[row_ord]
  abun_sel_ord <- abun_sel[row_ord, , drop = FALSE]

  col_ord <- order(abs(selected_df$z_main), decreasing = TRUE)
  abun_sel_ord <- abun_sel_ord[, col_ord, drop = FALSE]
  selected_df_ord <- selected_df[col_ord, , drop = FALSE]

  mat_plot <- log10(abun_sel_ord + 1e-6)
  direction_label <- ifelse(selected_df_ord$z_main < 0, "Higher in R", "Higher in NR")
  auxiliary_sum <- selected_df_ord$A + selected_df_ord$B + selected_df_ord$C

  annotation_col <- data.frame(
    Direction = factor(direction_label, levels = c("Higher in R", "Higher in NR")),
    "Auxiliary information" = auxiliary_sum,
    check.names = FALSE
  )
  rownames(annotation_col) <- colnames(mat_plot)

  annotation_row <- data.frame(
    Group = factor(as.character(outcome_ord), levels = c("NR", "R"))
  )
  rownames(annotation_row) <- rownames(mat_plot)

  gaps_row <- sum(outcome_ord == "NR")
  if (gaps_row <= 0 || gaps_row >= nrow(mat_plot)) {
    gaps_row <- NULL
  }

  n_features <- ncol(mat_plot)
  fontsize_col <- if (n_features <= 10) {
    12
  } else if (n_features <= 20) {
    10
  } else {
    8
  }
  cellwidth <- if (n_features <= 10) {
    32
  } else if (n_features <= 20) {
    24
  } else {
    18
  }
  left_padding_cm <- max(3.2, min(6.0, 0.15 * max(nchar(colnames(mat_plot)))))
  plot_width <- max(12, min(30, 6 + (cellwidth * n_features) / 72 + left_padding_cm / 2.54))
  plot_height <- max(8.5, min(11.5, 4.9 + 0.08 * nrow(mat_plot)))

  heatmap_obj <- pheatmap(
    mat = mat_plot,
    color = body_palette,
    breaks = body_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = TRUE,
    angle_col = 45,
    fontsize = 11,
    fontsize_col = fontsize_col,
    cellwidth = cellwidth,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    legend_breaks = legend_breaks,
    legend_labels = format(legend_breaks, trim = TRUE),
    gaps_row = gaps_row,
    main = title_text,
    silent = TRUE
  )

  gt <- heatmap_obj$gtable
  gt <- gtable::gtable_add_cols(gt, widths = grid::unit(left_padding_cm, "cm"), pos = 0)
  gt <- gtable::gtable_add_rows(gt, heights = grid::unit(0.45, "cm"), pos = 0)
  gt <- gtable::gtable_add_rows(gt, heights = grid::unit(0.7, "cm"), pos = -1)

  pdf_out <- file.path(output_dir, paste0("16s_LLM_family_", slug, "_weighted_heatmap.pdf"))
  png_out <- file.path(output_dir, paste0("16s_LLM_family_", slug, "_weighted_heatmap.png"))
  csv_out <- file.path(output_dir, paste0("16s_LLM_family_", slug, "_weighted_heatmap_features.csv"))

  write_padded_heatmap(gt, pdf_out, png_out, plot_width, plot_height)

  selected_df_ord %>%
    mutate(
      direction = direction_label,
      auxiliary_information = auxiliary_sum,
      display_rank = seq_len(n())
    ) %>%
    select(display_rank, taxon_key, z_main, p_main, A, B, C, auxiliary_information, direction) %>%
    write.csv(csv_out, row.names = FALSE)

  cat("Wrote:\n")
  cat(" - ", pdf_out, "\n", sep = "")
  cat(" - ", png_out, "\n", sep = "")
  cat(" - ", csv_out, "\n", sep = "")
}

write_heatmap_set(slice_head(ranked_df, n = 10), "top10", "Top 10 features by |z|")
write_heatmap_set(slice_head(ranked_df, n = 20), "top20", "Top 20 features by |z|")
write_heatmap_set(ranked_df, "all", "All features by |z|")
