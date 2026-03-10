library(devtools)
library(ape)
library(CATMicrobiome)
library(coin)
library(locfdr)
library(admix)
library(mvtnorm)
library(mixtools)
library(phyloseq)
library(devtools)
library(KScorrect)
library(mclust)
library(adaptMT)
library(AdaPTGMM)
library(splines2)
library(splines)
library(IHW)
library(FDRreg)

source("chai.R")
source("color_helper.R")

otuPath <- system.file("extdata","d1OTUtable.csv",
                       package = "CATMicrobiome")
otutable <- read.csv(otuPath,header=TRUE,row.names = 1)

taxonomyPath <- system.file("extdata","d1Taxonomy.csv",
                            package = "CATMicrobiome")
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

###################### Family level ps #######################
family <- tax_glom(ps, taxrank='family')


###################### Abundance table #######################
# family
cs <- colSums(otu_table(family))
cs[cs == 0] <- 1
family_abun <- sweep(otu_table(family), 2, cs, "/") * 100
t_family_abun <- t(family_abun)     # 43 x 76

########################## Outcome ###################################
# Check the metadata order matching to the abundance subject order
rownames(sample_data(family)) == rownames(t_family_abun)

outcome <- factor(sample_data(family)$BinOutcomes, levels = c("NR","R"))
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

# family's p-value and z-value
family_pz <- compute_pz_wilcox(t_family_abun, outcome)


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

family_pz_tax = make_main_df(t_abun = t_family_abun,  pz_df = family_pz,
                             ps_obj = family,  rank = "family")


#################### LLM results at family level ######################
library(dplyr)
library(readr)
library(auctestr)

class_aux_otu <- read_csv("./LLM_from_Yushu/d1Taxonomy_OTU_ABCD_class_and_below.csv")

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
        .cols  = c(A, B, C, D),
        .fns   = agg_fun,
        .names = "x{col}"
      ),
      .groups = "drop"
    ) %>%
    mutate(taxon_key = clean_taxon(.data[[level]]))

  aux_lvl
}

# # Family - Used class and below and aggregated at family/genus level
cla_ag_fam <- make_aux_tax_level(class_aux_otu, level = "family",
                                 agg_z_method = "stouffer")



############################ Fit Models ###############################
df_model <- family_pz_tax %>%
  dplyr::left_join(cla_ag_fam %>% dplyr::select(taxon_key, xA, xB, xC),
            by = "taxon_key")

z <- df_model$z_main
p <- df_model$p_main
X <- df_model[, c("xA", "xB", "xC")]

############################ Visualization ###############################
df <- data.frame(z, X)
df_long <- pivot_longer(df, cols = starts_with("x"), names_to = "variable", values_to = "value")
ggplot(df_long, aes(x = value, y = z)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~variable, scales = "free_x") +
  theme_minimal() +
  labs(title = "Relationship between z and predictors")

library(GGally)
p_llm_z <- ggpairs(df)

################## BH #######################
p_adjusted <- p.adjust(p, method = "BH")
bh_sel <- which(p_adjusted <= 0.05)

################## chai #######################
set.seed(123)
chai_family_llm <- chai(z, X, R = 100)

length(lFDRselect(chai_family_llm$lFDR, 0.05, 1)) # 14
length(lFDRselect(chai_family_llm$lFDR, 0.1, 1))  # 29

################## Adapt #######################
alphas <- seq(0.01, 0.10, by = 0.01)

# GLM
# ns
x <- data.frame(X)
formula_ns <- list(
  ~ ns(xA, df=2) + ns(xB, df=2) + ns(xC, df=2),
  ~ ns(xA, df=4) + ns(xB, df=4) + ns(xC, df=4),
  ~ ns(xA, df=6) + ns(xB, df=6) + ns(xC, df=6)
)

glm_ns <- adapt_glm(x, pvals = p, alphas=alphas,
                    pi_formulas = formula_ns, mu_formulas = formula_ns)

# GMM - p
# ns
gmm_ns_p <- adapt_gmm(x, pvals = p,  alphas=alphas,
                      beta_formulas = formula_ns)
# alpha = 0.1: FDPhat 0.075, Number of Rej. 2
# alpha = 0.05: FDPhat 0.05, Number of Rej. 1

# GMM - Z
gmm_ns_z <- adapt_gmm(x, z = z, alphas=alphas,
                      beta_formulas = formula_ns, testing = "two_sided")
# alpha = 0.1: FDPhat 0.075, Number of Rej. 0
# alpha = 0.05: FDPhat 0.05, Number of Rej. 0


################## OrderShapeEM #######################
require(OrderShapeEM)
ordershape_x1 <- OrderShapeEM(pvals = p, order.var = X[,1],
                              OrderShapeEM.control(trace = TRUE))
sum(ordershape_x1$fdr <= 0.05)
# 3

ordershape_x2 <- OrderShapeEM(pvals = p, order.var = X[,2],
                              OrderShapeEM.control(trace = TRUE))
sum(ordershape_x2$fdr <= 0.05)
# 28

ordershape_x3 <- OrderShapeEM(pvals = p, order.var = X[,3],
                              OrderShapeEM.control(trace = TRUE))
sum(ordershape_x3$fdr <= 0.05)
# 4


######################## FDRreg ############################
set.seed(123)
fdr_theo <- FDRreg(z, as.matrix(X), nulltype = 'theoretical')
length(which(fdr_theo$FDR <= 0.05))
# 3

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
  length(lFDRselect(chai_family_llm$lFDR, q, 1))
})

# adapt_glm
adapt_glm_ns_rej <- glm_ns$nrejs[1:10]

# adapt_gmm_p
adapt_gmm_p_ns_rej <- gmm_ns_p$nrejs[1:10]
# adapt_gmm_p_bs_rej <- gmm_bs_p$nrejs[1:10]

# adapt_gmm_z
adapt_gmm_z_ns_rej <- gmm_ns_z$nrejs[1:10]

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
rdrreg_rejs <- sapply(q_levels, function(q) {
  length(which(fdr_theo$FDR <= q))
})

# IHW
# IHW
ihw_1_rejs <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p, covariates = X[,1], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_2_rejs <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p, covariates = X[,2], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_3_rejs <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p, covariates = X[,3], alpha = q, nbins = 5)
  rejections(ihw)
})

# BH
p_adjusted <- p.adjust(p, method = "BH")

bh_rejs <- sapply(q_levels, function(q) {
  length(which(p_adjusted <= q))
})

# combine into data frame
LLM_fam_classandbelow <- data.frame(
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




###########################################################
# Plots
library(ggplot2)
library(tidyr)
library(ggrepel)
library(viridisLite)

# into long format
LLM_fam_classandbelow_long <- LLM_fam_classandbelow %>%
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


plot_rejections_vs_q <- function(long_table, title = "Number of Discoveries vs q",
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
    labs(title = title, x = "q level", y = "Number of Rejections") +
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

p_LLM <- plot_rejections_vs_q(LLM_fam_classandbelow_long,
                                  title = "Number of Discoveries vs q, with LLM",
                                  x_breaks = sort(unique(LLM_fam_classandbelow_long$q)))

# pdf("./plots/Gopalakrishnan/16s_wilcoxz_LLM_fam_classandbelow.pdf", width = 10, height = 5)
# plot(p_LLM)
# dev.off()



# Visualization related to chai results
chai_sel <- lFDRselect(chai_family_llm$lFDR, 0.05, 1)
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

ggplot(df_long, aes(x = value, y = z)) +
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
  labs(title = "Relationship between z and predictors")

ggplot(df, aes(x = xA, y = xB, color = grp)) +
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
    title = "LLM Space Map",
    subtitle = "Points colored by findings",
    x = "LLM 1",
    y = "LLM 2",
    color = "Group"
  ) +
  theme(legend.position = "right")

# Explore the selections
library(ape); library(phyloseq)

# inputs
tr <- phy_tree(family)                 # phylo
sel_otus <- df_model$otu_rep[chai_sel]   # chai selected OTU IDs
# keep only OTUs that exist in the tree (safe)
sel_otus <- intersect(sel_otus, tr$tip.label)

bh_otus <- df_model$otu_rep[bh_sel]   # BH selected OTU IDs


# prune tree to selected OTUs
tr_sel <- keep.tip(tr, sel_otus)

# Brownian-motion phylogenetic covariance, then correlation
V <- vcv.phylo(tr_sel)     # covariance implied by tree branch lengths
R_phylo <- cov2cor(V)      # correlation matrix

# Rename OTU to taxa
tax <- as.data.frame(tax_table(family))
tax$OTU <- rownames(tax)
lab_df <- tax %>%
  filter(OTU %in% sel_otus) %>%
  mutate(
    raw = .data[["family"]],
    clean = clean_taxon(raw),
    label = ifelse(is.na(clean) | clean == "", OTU, paste0(clean, " (", OTU, ")"))
  )

lbl <- setNames(make.unique(lab_df$label), lab_df$OTU)

rownames(R_phylo) <- lbl[rownames(R_phylo)]
colnames(R_phylo) <- lbl[colnames(R_phylo)]
R_phylo[1:5, 1:5]

library(pheatmap)
pheatmap(
  R_phylo,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "Phylogenetic correlation among selected taxa"
)


# Tree
tr_all <- phy_tree(family)

# map labels for ALL tips (not only selected)
tax_all <- as.data.frame(tax_table(family))
tax_all$OTU <- rownames(tax_all)

rank_to_show <- "family"

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
plot(tr_all, show.tip.label = FALSE, cex = 0.7, tip.color = tip_cols,
     main = "Phylogenetic tree with selected families by chai using LLM")

# add cleaned labels on top (so you see taxa names)
tiplabels(text = unname(lbl_all[tr_all$tip.label]) , tip = 1:length(tr_all$tip.label),
          frame = "none", adj = 0, cex = 0.7, col = tip_cols)


# Circular tree
library(dplyr)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

rank_to_show <- "family"

anno <- tax_all %>%
  # tibble::rownames_to_column("OTU") %>%
  mutate(
    # 清洗要显示的标签（family）
    raw = .data[[rank_to_show]],
    clean = clean_taxon(raw),
    tip_label = ifelse(is.na(clean) | clean == "", OTU, clean),

    # 清洗 phylum（用于色块）
    phylum_clean = clean_taxon(phylum),

    # tip 文字颜色（红/黑）
    tip_color = dplyr::case_when(
      OTU %in% sel_otus & OTU %in% bh_otus ~ "both",
      OTU %in% sel_otus & !(OTU %in% bh_otus) ~ "chai_only",
      !(OTU %in% sel_otus) & OTU %in% bh_otus ~ "BH_only",
      TRUE ~ "neither"
    )  ) %>%
  dplyr::select(OTU, tip_label, phylum_clean, tip_color) %>%
  distinct(OTU, .keep_all = TRUE)     # 确保顺序与树 tip 一致（ggtree 依赖 label 对齐）

anno$tip_color <- factor(
  anno$tip_color,
  levels = c("chai_only", "BH_only", "both", "neither")
)

tr <- ape::as.phylo(tr_all)
p <- ggtree(tr, layout = "circular", size = 0.5)

p1 <- p %<+% anno

# 外圈 phylum 色块（每个 tip 一个小方块）
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

# 先关掉 tiplab 的图例（避免它贡献 legend key）
p3 <- p2 +
  geom_tiplab(
    aes(label = tip_label, color = tip_color),
    size = 2,
    offset = 0.02,          # 文字离色块距离：想更远就调大
    show.legend = FALSE
  ) +
  scale_color_manual(
    name   = "Tip label group",
    values = c(chai_only="red", BH_only="blue", both="green", neither="black"),
    breaks = c("chai_only", "BH_only", "both", "neither"),
    labels = c("Only in chai", "Only in BH", "In both", "Neither"),
    drop   = FALSE
  ) +
  labs(title = "Phylogenetic tree with selected families by chai on wilcoxon z with LLM") +
  theme(legend.position = "right")

# 2) 加一个“只为图例服务”的点层（不会画在图上，但会生成完美图例）
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
#
# pdf("./plots/Gopalakrishnan/16s_wilcoxz_LLM_phylo_tree.pdf", width = 8, height = 6.5)
# plot(p4)
# dev.off()


################# Check which features have higher influence ###########################
full <- lab_all %>%
  left_join(cla_ag_fam, by = c("label" = "family"))

top_sig_A <- full %>%
  dplyr::arrange(desc(abs(xA))) %>%
  head(20)
print(top_sig_A)

top_sig_B <- full %>%
  dplyr::arrange(desc(abs(xB))) %>%
  head(20)
print(top_sig_B)

top_sig_C <- full %>%
  dplyr::arrange(desc(abs(xC))) %>%
  head(20)
print(top_sig_C)

