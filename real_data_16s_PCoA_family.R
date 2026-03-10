##################################### family ########################################
library(devtools)
library(ape)
library(CATMicrobiome)
library(phyloseq)
library(coin)

# load("C:/Users/zwang26/OneDrive - UTHealth Houston/conditionalGaussian/chai_env/Gopalakrishnan_16s_tables.RData")


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

###################### Family level #######################
family <- tax_glom(ps, taxrank='family')
# genus_vcv <- vcv(phy_tree(genus))
family_D <- cophenetic(phy_tree(family))
# 76 x 76

###################### PCoA ##########################
family_pcoa <- cmdscale(family_D, k = 2)
table(rownames(tax_table(family)) == rownames(family_pcoa))


# Convert to relative abundance
cs <- colSums(otu_table(family))
cs[cs == 0] <- 1
family_abun <- sweep(otu_table(family), 2, cs, "/") * 100    # family x subjects

# Transpose for modeling (subjects × species)
t_abun <- t(family_abun)     # 43 x 76


########################## Outcome ###################################
# Check the metadata order matching to the abundance subject order
rownames(sample_data(family)) == rownames(t_abun)

outcome <- factor(sample_data(family)$BinOutcomes, levels = c("NR","R"))
table(outcome)


##################### p-value/z-value ###########################
p_z_values <- sapply(colnames(t_abun), function(feature) {
  features <- t_abun[, feature, drop = FALSE]
  wilcox <- wilcox_test(features[,1] ~ outcome)

  p_value <- pvalue(wilcox)
  z_value <- statistic(wilcox)

  return(c(p_value = p_value, z_value = z_value))
}, simplify = "matrix")
p_z_values <- t(p_z_values)

p_value <- p_z_values[,1]
z_value <- p_z_values[,2]

##################### BH ###########################
p_adjusted <- p.adjust(p_value, method = "BH")
bh <- which(p_adjusted < 0.05)
length(bh) # 1





########################### Plot of z and X ##############################
df <- data.frame(z = z_value, x = family_pcoa)
# z
ggplot(df, aes(x=z)) +
  geom_histogram(aes(y=..density..), binwidth=0.25, color="black", fill="lightblue") +
  geom_density() +
  stat_function(fun=dnorm, args=list(mean=0, sd=1), linetype="dashed") +
  labs(title="Histogram of z-value with N(0,1)", x="z", y="Density") +
  theme_minimal()

# x - PCoA 1st
ggplot(df, aes(x = x.1)) +
  geom_histogram(aes(y = ..density..), bins = 40, color = "black", fill = "lightblue") +
  geom_density() +
  labs(
    title = paste0("Histogram of X - PCoA 1st dimension"),
    x = "X", y = "Density"
  ) +
  theme_minimal()

# x - PCoA 2nd
ggplot(df, aes(x = x.2)) +
  geom_histogram(aes(y = ..density..), bins = 40, color = "black", fill = "lightblue") +
  geom_density() +
  labs(
    title = paste0("Histogram of X - PCoA 2nd dimension"),
    x = "X", y = "Density"
  ) +
  theme_minimal()

# z vs. x
# ggplot(df, aes(x = z, y = x.1)) +
#   geom_point(alpha = 0.6) +
#   labs(x = "z-value", y = "Canonical Correlation") +
#   theme_minimal()

# p
ggplot(data.frame(p_value), aes(p_value)) +
  geom_histogram(bins=40, color="black", fill="grey80") +
  labs(title="Histogram of two-sided p-values", x="p-value", y="Count") +
  theme_minimal()

############################ Fit model ############################
library(locfdr)
library(admix)
library(mvtnorm)
library(mixtools)
library(devtools)
library(KScorrect)
library(mclust)
library(adaptMT)
library(AdaPTGMM)
library(splines2)
library(splines)
library(IHW)
library(FDRreg)

setwd("C:/Users/zwang26/OneDrive - UTHealth Houston/conditionalGaussian")
# source("ConditionalGaussianlFDR_defined_functions_playing.R")
source("chai.R")
source("color_helper.R")


################## chai #######################
set.seed(123)
chai_family_pcoa <- chai(z_value, family_pcoa, M = 100)

length(lFDRselect(chai_family_pcoa$lFDR, 0.05, 1)) # 15
length(lFDRselect(chai_family_pcoa$lFDR, 0.1, 1))  # 31


# Visualization
chai_pcoa_sel <- lFDRselect(chai_family_pcoa$lFDR, 0.05, 1)
df <- data.frame(z_value, family_pcoa) %>%
  dplyr::mutate(
    grp = dplyr::case_when(
      dplyr::row_number() %in% bh & dplyr::row_number() %in% chai_pcoa_sel ~ "both",
      dplyr::row_number() %in% bh ~ "BH_only",
      dplyr::row_number() %in% chai_pcoa_sel ~ "chai_only",
      TRUE ~ "neither"
    )
  )

ggplot(df, aes(x = family_pcoa[,1], y = family_pcoa[,2], color = grp)) +
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
    title = "PCoA of Phylogenatic Tree Structure Space Map",
    subtitle = "Points colored by findings",
    x = "PCoA Dimension 1",
    y = "PCoA Dimension 2",
    color = "Group"
  ) +
  theme(legend.position = "right")


# map labels for ALL tips (not only selected)
tax_all <- as.data.frame(tax_table(family))
tax_all$OTU <- rownames(tax_all)

tax_all[bh,]


################## Adapt #######################
alphas <- seq(0.01, 0.10, by = 0.01)

# GLM
# ns
x <- data.frame(family_pcoa)
formula_ns <- list("~ns(X1, df = 2) + ns(X2, df = 2)",
                   "~ns(X1, df = 4) + ns(X2, df = 4)",
                   "~ns(X1, df = 6) + ns(X2, df = 6)")
glm_ns <- adapt_glm(x, pvals = p_value, alphas=alphas,
                    pi_formulas = formula_ns, mu_formulas = formula_ns)
# alpha = 0.1: FDPhat 0.0833, Number of Rej. 12
# alpha = 0.09: FDPhat 0.0833, Number of Rej. 12

# # ispline
# x <- data.frame(family_pcoa)
# formula_is <- list("~iSpline(X1, df = 4) + iSpline(X2, df = 4)",
#                    "~iSpline(X1, df = 6) + iSpline(X2, df = 6)",
#                    "~iSpline(X1, df = 8) + iSpline(X2, df = 8)")
# glm_ispline <- adapt_glm(x, pvals = p_value,  alphas=alphas,
#                          pi_formulas = formula_is, mu_formulas = formula_is)
# # alpha = 0.1: FDPhat 0.0833, Number of Rej. 12
# # alpha = 0.09: FDPhat 0.0833, Number of Rej. 12
#
# # bs
# formula_bs <- list("~bs(X1, df = 4) + bs(X2, df = 4)",
#                    "~bs(X1, df = 6) + bs(X2, df = 6)",
#                    "~bs(X1, df = 8) + bs(X2, df = 8)")
# glm_bs <- adapt_glm(x, pvals = p_value,  alphas=alphas,
#                     pi_formulas = formula_bs, mu_formulas = formula_bs)
# # 0

# GMM - p
# ns
gmm_ns_p <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
                      beta_formulas = formula_ns)
# alpha = 0.1: FDPhat 0.075, Number of Rej. 2
# alpha = 0.05: FDPhat 0.05, Number of Rej. 1


# # ispline
# gmm_ispline_p <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
#                            beta_formulas = formula_is)
# # alpha = 0.1: FDPhat 0.075, Number of Rej. 2
# # alpha = 0.06: FDPhat 0.05, Number of Rej. 1
#
# # bs
# gmm_bs_p <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
#                       beta_formulas = formula_bs)
#
# # alpha = 0.1: FDPhat 0.075, Number of Rej. 2
# # alpha = 0.05: FDPhat 0.05, Number of Rej. 1

# GMM - Z
gmm_ns_z <- adapt_gmm(x, z = z_value, alphas=alphas,
                      beta_formulas = formula_ns, testing = "two_sided")
# alpha = 0.1: FDPhat 0.075, Number of Rej. 2
# alpha = 0.05: FDPhat 0.05, Number of Rej. 2

# gmm_bs_z <- adapt_gmm(x, z = z_value, alphas=alphas,
#                       beta_formulas = formula_bs, testing = "two_sided")
# # alpha = 0.1: FDPhat 0.4429, Number of Rej. 0
# # alpha = 0.05: FDPhat 0.4429, Number of Rej. 0

################## OrderShapeEM #######################
require(OrderShapeEM)
ordershape_x1 <- OrderShapeEM(pvals = p_value, order.var = family_pcoa[,1],
                              OrderShapeEM.control(trace = TRUE))
sum(ordershape_x1$fdr <= 0.05)
# 2

ordershape_x2 <- OrderShapeEM(pvals = p_value, order.var = family_pcoa[,2],
                              OrderShapeEM.control(trace = TRUE))
sum(ordershape_x2$fdr <= 0.05)
# 11

######################## FDRreg ############################
set.seed(123)
fdr_theo <- FDRreg(z_value, family_pcoa, nulltype = 'theoretical')
length(which(fdr_theo$FDR <= 0.05))
# 4

fdr_emp <- FDRreg(z_value, family_pcoa, nulltype = 'empirical')
length(which(fdr_emp$FDR <= 0.05))
# 14

######################## IHW ############################
set.seed(123)
ihw_x1_01 <- ihw(pvalues = p_value, covariates = family_pcoa[,1], alpha = 0.01, nbins = 5)
ihw_x1_02 <- ihw(pvalues = p_value, covariates = family_pcoa[,1], alpha = 0.02, nbins = 5)
ihw_x1_03 <- ihw(pvalues = p_value, covariates = family_pcoa[,1], alpha = 0.03, nbins = 5)
ihw_x1_04 <- ihw(pvalues = p_value, covariates = family_pcoa[,1], alpha = 0.04, nbins = 5)
ihw_x1_05 <- ihw(pvalues = p_value, covariates = family_pcoa[,1], alpha = 0.05, nbins = 5)
ihw_x1_06 <- ihw(pvalues = p_value, covariates = family_pcoa[,1], alpha = 0.06, nbins = 5)
ihw_x1_07 <- ihw(pvalues = p_value, covariates = family_pcoa[,1], alpha = 0.07, nbins = 5)
ihw_x1_08 <- ihw(pvalues = p_value, covariates = family_pcoa[,1], alpha = 0.08, nbins = 5)
ihw_x1_09 <- ihw(pvalues = p_value, covariates = family_pcoa[,1], alpha = 0.09, nbins = 5)
ihw_x1_1 <- ihw(pvalues = p_value, covariates = family_pcoa[,1], alpha = 0.1, nbins = 5)

ihw_x2_01 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.01, nbins = 5)
ihw_x2_02 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.02, nbins = 5)
ihw_x2_03 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.03, nbins = 5)
ihw_x2_04 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.04, nbins = 5)
ihw_x2_05 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.05, nbins = 5)
ihw_x2_06 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.06, nbins = 5)
ihw_x2_07 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.07, nbins = 5)
ihw_x2_08 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.08, nbins = 5)
ihw_x2_09 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.09, nbins = 5)
ihw_x2_1 <- ihw(pvalues = p_value, covariates = family_pcoa[,2], alpha = 0.1, nbins = 5)

rejections(ihw_x1_05)  # 1   1
rejections(ihw_x1_06)  # 1   1
rejections(ihw_x1_07)  # 1   1
rejections(ihw_x1_08)  # 0   1
rejections(ihw_x1_09)  # 0   1
rejections(ihw_x1_1)   # 0   2

###########################################################
# Summarize into a different q level table
q_levels <- alphas

# chai
chai_rejs <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_family_pcoa$lFDR, q, 1))
})

# adapt_glm
adapt_glm_ns_rej <- glm_ns$nrejs[1:10]
# adapt_glm_bs_rej <- glm_bs$nrejs[1:10]

# adapt_gmm_p
adapt_gmm_p_ns_rej <- gmm_ns_p$nrejs[1:10]
# adapt_gmm_p_bs_rej <- gmm_bs_p$nrejs[1:10]

# adapt_gmm_z
adapt_gmm_z_ns_rej <- gmm_ns_z$nrejs[1:10]
# adapt_gmm_z_bs_rej <- gmm_bs_z$nrejs[1:10]

# OrderShapeEM
ordershapeem_x1_rejs <- sapply(q_levels, function(q) {
  sum(ordershape_x1$fdr <= q)
})

ordershapeem_x2_rejs <- sapply(q_levels, function(q) {
  sum(ordershape_x2$fdr <= q)
})

# FDRreg
rdrreg_rejs <- sapply(q_levels, function(q) {
  length(which(fdr_theo$FDR <= q))
})

# IHW
ihw_x1_rejs <- c(rejections(ihw_x1_01), rejections(ihw_x1_02), rejections(ihw_x1_03),
                 rejections(ihw_x1_04),rejections(ihw_x1_05), rejections(ihw_x1_06),
                 rejections(ihw_x1_07), rejections(ihw_x1_08),
                 rejections(ihw_x1_09), rejections(ihw_x1_1))
# [1] 0 1 1 1 1 1 1 0 0 0

ihw_x2_rejs <- c(rejections(ihw_x2_01), rejections(ihw_x2_02), rejections(ihw_x2_03),
                 rejections(ihw_x2_04), rejections(ihw_x2_05), rejections(ihw_x2_06),
                 rejections(ihw_x2_07), rejections(ihw_x2_08),
                 rejections(ihw_x2_09), rejections(ihw_x2_1))
# [1] 0 1 1 1 1 1 1 1 1 2

# BH
bh_rejs <- sapply(q_levels, function(q) {
  length(which(p_adjusted <= q))
})
# [1] 0 1 1 1 1 1 1 1 1 2



# combine into data frame
gopa_sum_16s_family <- data.frame(
  q = q_levels,
  chai = chai_rejs,
  adapt_glm = adapt_glm_ns_rej,
  # adapt_glm_bs = adapt_glm_bs_rej,
  adapt_gmm_p = adapt_gmm_p_ns_rej,
  # adapt_gmm_p_bs = adapt_gmm_p_bs_rej,
  adapt_gmm_z = adapt_gmm_z_ns_rej,
  # adapt_gmm_z_bs = adapt_gmm_z_bs_rej,
  OrderShapeEM_x1 = ordershapeem_x1_rejs,
  OrderShapeEM_x2 = ordershapeem_x2_rejs,
  FDRreg = rdrreg_rejs,
  IHW_x1 = ihw_x1_rejs,
  IHW_x2 = ihw_x2_rejs,
  BH = bh_rejs
)



###########################################################
# Plots
library(ggplot2)
library(tidyr)
library(ggrepel)
library(viridisLite)

# into long format
gopa_16s_long <- gopa_sum_16s_family %>%
  pivot_longer(
    cols = -q,
    names_to = "Method",
    values_to = "Rejections"
  ) %>%
  mutate(Method = as.character(Method))


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

p_gopa_16s_pcoa <- plot_rejections_vs_q(gopa_16s_long,
                                   title = "Number of Discoveries vs q, on Wilcoxon z-statistics",
                                   x_breaks = sort(unique(gopa_16s_long$q)))

pdf("./plots/Gopalakrishnan/16s_wilcoxz_pcoa_family_updated.pdf", width = 8, height = 4)
plot(p_gopa_16s_pcoa)
dev.off()


# # into long format
# eraw_long_16s_family <- eraw_sum_16s_family %>%
#   pivot_longer(
#     cols = -q,
#     names_to = "Method",
#     values_to = "Rejections"
#   )
#
# eraw_chai <- filter(eraw_long_16s_family, Method == "chai")
# eraw_others <- filter(eraw_long_16s_family, Method != "chai")
#
# # install.packages("ggrepel")  # once
# library(ggrepel)
# library(viridisLite)
#
# # color map: viridis for others, red for chai
# methods <- unique(eraw_long_16s_family$Method)
# others  <- setdiff(methods, "chai")
# pal     <- viridis(length(others), option = "plasma")
# method_cols <- c(setNames(pal, others), "chai" = "red")
#
# # last point (max q) for each method
# endpoints <- eraw_long_16s_family |>
#   dplyr::arrange(Method, q) |>
#   dplyr::group_by(Method) |>
#   dplyr::slice_tail(n = 1) |>
#   dplyr::ungroup()
#
# # plot with end-of-line labels
# p <- ggplot() +
#   geom_line(data = eraw_others, aes(q, Rejections, color = Method), size = 0.9) +
#   geom_point(data = eraw_others, aes(q, Rejections, color = Method), size = 1.8) +
#   geom_line(data = eraw_chai,   aes(q, Rejections, color = Method), size = 1.5) +
#   geom_point(data = eraw_chai,  aes(q, Rejections, color = Method), size = 2.5) +
#   geom_text_repel(
#     data = endpoints,
#     aes(q, Rejections, label = Method, color = Method),
#     nudge_x = diff(range(eraw_long_16s_family$q)) * 0.03,  # push labels a bit to the right
#     hjust = 0, direction = "y",
#     segment.size = 0.2, box.padding = 0.1, point.padding = 0.1,
#     show.legend = FALSE
#   ) +
#   scale_color_manual(values = method_cols) +
#   labs(title = "Number of Discoveries vs q", x = "q level", y = "Number of Rejections") +
#   coord_cartesian(clip = "off") +
#   expand_limits(x = max(eraw_long_16s_family$q) + 0.01) +         # room for labels
#   theme_minimal() +
#   theme(
#     plot.margin = margin(5.5, 60, 5.5, 5.5),           # extra right margin
#     legend.position = "none"                            # labels replace legend; set to "bottom" if you want both
#   )










########################################################################################################
############################## DESeq2 ######################################
library(DESeq2)
set.seed(123)
cts <- as(otu_table(family), "matrix")

md  <- as(sample_data(family), "data.frame")
md$BinOutcomes <- factor(md$BinOutcomes, levels = c("NR", "R"))

dds <- DESeqDataSetFromMatrix(countData = cts, colData = md,
                              design = ~ BinOutcomes)
dds_res <- DESeq(dds)
res <- results(dds_res, alpha = 0.05)
summary(res)
which(res$padj <= 0.05)
# [1]  9 28 48 54 63 65 76

# res_01 <- results(dds_res, alpha = 0.01)  # 6
# res_02 <- results(dds_res, alpha = 0.02)  # 6
# res_03 <- results(dds_res, alpha = 0.03)  # 7
# res_04 <- results(dds_res, alpha = 0.04)  # 7
res_05 <- results(dds_res, alpha = 0.05)  # 7
# res_06 <- results(dds_res, alpha = 0.06)  # 7
# res_07 <- results(dds_res, alpha = 0.07)  # 8
# res_08 <- results(dds_res, alpha = 0.08)  # 8
# res_09 <- results(dds_res, alpha = 0.09)  # 8
# res_1 <- results(dds_res, alpha = 0.1)    # 8

DESeq <- c(6,6,7,7,7,7,8,8,8,8)

# res_1 <- results(dds_res, alpha = 0.1)
# summary(res_1)
# which(res_1$padj <= 0.1)
# # [1]  9 27 28 48 54 63 65 76


DES_stat <- res_05$stat
DES_p <- res_05$pvalue

########################### BH ###################################
DES_p_adjusted <- p.adjust(DES_p, method = "BH")
DES_bh <- which(DES_p_adjusted <= 0.05)
length(DES_bh) # 7
# [1] "OTU_175" "OTU_83"  "OTU_16"  "OTU_95"  "OTU_86"  "OTU_7"   "OTU_3"

########################### Plot of z and X ##############################
df <- data.frame(z = DES_stat, x = family_pcoa)
# z
ggplot(df, aes(x=z)) +
  geom_histogram(aes(y=..density..), binwidth=0.25, color="black", fill="lightblue") +
  geom_density() +
  stat_function(fun=dnorm, args=list(mean=0, sd=1), linetype="dashed") +
  labs(title="Histogram of DESeq2 Wald statistic with N(0,1)", x="z", y="Density") +
  theme_minimal()


# p
ggplot(data.frame(DES_p), aes(DES_p)) +
  geom_histogram(bins=40, color="black", fill="grey80") +
  labs(title="Histogram of two-sided p-values", x="p-value", y="Count") +
  theme_minimal()

############################ Fit model ############################
################## chai #######################
set.seed(123)
chai_DES_pcoa <- chai(DES_stat, family_pcoa, R = 100)

length(lFDRselect(chai_DES_pcoa$lFDR, 0.05, 1)) # 14
length(lFDRselect(chai_DES_pcoa$lFDR, 0.1, 1))  # 23

# Quickly check which family taxa were been picked
chai_found <- tax_table(family)[rownames(tax_table(family)) %in% rownames(family_pcoa)[lFDRselect(chai_DES_pcoa$lFDR, 0.05, 1)], ]
# write.csv(chai_found, "C:/Users/zwang26/OneDrive - UTHealth Houston/conditionalGaussian/plots/Gopalakrishnan/chai_found_DES_q05_family.csv", row.names = TRUE)

# Visualization
chai_DES_sel <- lFDRselect(chai_DES_pcoa$lFDR, 0.05, 1)
df_DES <- data.frame(DES_stat, family_pcoa) %>%
  dplyr::mutate(
    grp = dplyr::case_when(
      dplyr::row_number() %in% DES_bh & dplyr::row_number() %in% chai_DES_sel ~ "both",
      dplyr::row_number() %in% DES_bh ~ "BH_only",
      dplyr::row_number() %in% chai_DES_sel ~ "chai_only",
      TRUE ~ "neither"
    )
  )

ggplot(df_DES, aes(x = family_pcoa[,1], y = family_pcoa[,2], color = grp)) +
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
    title = "PCoA of Phylogenatic Tree Structure Space Map",
    subtitle = "Points colored by findings",
    x = "PCoA Dimension 1",
    y = "PCoA Dimension 2",
    color = "Group"
  ) +
  theme(legend.position = "right")


################## Adapt #######################
alphas <- seq(0.01, 0.10, by = 0.01)

# GLM
# ns
x <- data.frame(family_pcoa)

glm_DES_ns <- adapt_glm(x, pvals = DES_p, alphas=alphas,
                        pi_formulas = formula_ns, mu_formulas = formula_ns)
# 0

# # ispline
# x <- data.frame(family_pcoa)
#
# glm_DES_ispline <- adapt_glm(x, pvals = DES_p,  alphas=alphas,
#                              pi_formulas = formula_is, mu_formulas = formula_is)
# # 0

# bs

# glm_DES_bs <- adapt_glm(x, pvals = DES_p,  alphas=alphas,
#                         pi_formulas = formula_bs, mu_formulas = formula_bs)
# # 0

# GMM - p
# ns
gmm_DES_ns_p <- adapt_gmm(x, pvals = DES_p,  alphas=alphas,
                          beta_formulas = formula_ns)
# alpha = 0.1: FDPhat 0.0944, Number of Rej. 9
# alpha = 0.05: FDPhat 0.05, Number of Rej. 9


# # ispline
# gmm_DES_ispline_p <- adapt_gmm(x, pvals = DES_p,  alphas=alphas,
#                                beta_formulas = formula_is)
# # alpha = 0.1: FDPhat 0.0938, Number of Rej. 8
# # alpha = 0.05: FDPhat 0.05, Number of Rej. 7

# # bs
# gmm_DES_bs_p <- adapt_gmm(x, pvals = DES_p,  alphas=alphas,
#                           beta_formulas = formula_bs)
# # alpha = 0.1: FDPhat 0.0938, Number of Rej. 8
# # alpha = 0.05: FDPhat 0.05, Number of Rej. 7


# GMM - Z
gmm_DES_ns_z <- adapt_gmm(x, z = DES_stat, alphas=alphas,
                          beta_formulas = formula_ns, testing = "two_sided")
# alpha = 0.1: FDPhat 0.0938, Number of Rej. 8
# alpha = 0.05: FDPhat 0.05, Number of Rej. 7

# gmm_DES_bs_z <- adapt_gmm(x, z = DES_stat, alphas=alphas,
#                           beta_formulas = formula_bs, testing = "two_sided")
# # alpha = 0.1: FDPhat 0.0938, Number of Rej. 8
# # alpha = 0.05: FDPhat 0.05, Number of Rej. 7

################## OrderShapeEM #######################

require(OrderShapeEM)
ordershape_DES_x1 <- OrderShapeEM(pvals = DES_p, order.var = family_pcoa[,1],
                                  OrderShapeEM.control(trace = TRUE))
sum(ordershape_DES_x1$fdr <= 0.05)
# 0

ordershape_DES_x2 <- OrderShapeEM(pvals = DES_p, order.var = family_pcoa[,2],
                                  OrderShapeEM.control(trace = TRUE))
sum(ordershape_DES_x2$fdr <= 0.05)
# 0

######################## FDRreg ############################
set.seed(123)
fdr_DES_theo <- FDRreg(DES_stat, family_pcoa, nulltype = 'theoretical')
length(which(fdr_DES_theo$FDR <= 0.05))
# 9

fdr_DES_emp <- FDRreg(DES_stat, family_pcoa, nulltype = 'empirical')
length(which(fdr_DES_emp$FDR <= 0.05))
# 9



######################## IHW ############################
ihw_DES_x1_01 <- ihw(pvalues = DES_p, covariates = family_pcoa[,1], alpha = 0.01, nbins = 5)
ihw_DES_x1_02 <- ihw(pvalues = DES_p, covariates = family_pcoa[,1], alpha = 0.02, nbins = 5)
ihw_DES_x1_03 <- ihw(pvalues = DES_p, covariates = family_pcoa[,1], alpha = 0.03, nbins = 5)
ihw_DES_x1_04 <- ihw(pvalues = DES_p, covariates = family_pcoa[,1], alpha = 0.04, nbins = 5)
ihw_DES_x1_05 <- ihw(pvalues = DES_p, covariates = family_pcoa[,1], alpha = 0.05, nbins = 5)
ihw_DES_x1_06 <- ihw(pvalues = DES_p, covariates = family_pcoa[,1], alpha = 0.06, nbins = 5)
ihw_DES_x1_07 <- ihw(pvalues = DES_p, covariates = family_pcoa[,1], alpha = 0.07, nbins = 5)
ihw_DES_x1_08 <- ihw(pvalues = DES_p, covariates = family_pcoa[,1], alpha = 0.08, nbins = 5)
ihw_DES_x1_09 <- ihw(pvalues = DES_p, covariates = family_pcoa[,1], alpha = 0.09, nbins = 5)
ihw_DES_x1_1 <- ihw(pvalues = DES_p, covariates = family_pcoa[,1], alpha = 0.1, nbins = 5)

ihw_DES_x2_01 <- ihw(pvalues = DES_p, covariates = family_pcoa[,2], alpha = 0.01, nbins = 5)
ihw_DES_x2_02 <- ihw(pvalues = DES_p, covariates = family_pcoa[,2], alpha = 0.02, nbins = 5)
ihw_DES_x2_03 <- ihw(pvalues = DES_p, covariates = family_pcoa[,2], alpha = 0.03, nbins = 5)
ihw_DES_x2_04 <- ihw(pvalues = DES_p, covariates = family_pcoa[,2], alpha = 0.04, nbins = 5)
ihw_DES_x2_05 <- ihw(pvalues = DES_p, covariates = family_pcoa[,2], alpha = 0.05, nbins = 5)
ihw_DES_x2_06 <- ihw(pvalues = DES_p, covariates = family_pcoa[,2], alpha = 0.06, nbins = 5)
ihw_DES_x2_07 <- ihw(pvalues = DES_p, covariates = family_pcoa[,2], alpha = 0.07, nbins = 5)
ihw_DES_x2_08 <- ihw(pvalues = DES_p, covariates = family_pcoa[,2], alpha = 0.08, nbins = 5)
ihw_DES_x2_09 <- ihw(pvalues = DES_p, covariates = family_pcoa[,2], alpha = 0.09, nbins = 5)
ihw_DES_x2_1 <- ihw(pvalues = DES_p, covariates = family_pcoa[,2], alpha = 0.1, nbins = 5)

rejections(ihw_DES_x1_05)  # 7   7
rejections(ihw_DES_x1_06)  # 7   7
rejections(ihw_DES_x1_07)  # 8   8
rejections(ihw_DES_x1_08)  # 8   8
rejections(ihw_DES_x1_09)  # 8   8
rejections(ihw_DES_x1_1)   # 8   8

###########################################################
# Summarize into a different q level table
q_levels <- alphas

# chai
chai_DES_rejs <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_DES_pcoa$lFDR, q, 1))
})

# adapt_glm
adapt_glm_DES_ns_rej <- glm_DES_ns$nrejs[1:10]
# adapt_glm_DES_bs_rej <- glm_DES_bs$nrejs[1:10]

# adapt_gmm_p
adapt_gmm_DES_p_ns_rej <- gmm_DES_ns_p$nrejs[1:10]
# adapt_gmm_DES_p_bs_rej <- gmm_DES_bs_p$nrejs[1:10]

# adapt_gmm_z
adapt_gmm_DES_z_ns_rej <- gmm_DES_ns_z$nrejs[1:10]
# adapt_gmm_DES_z_bs_rej <- gmm_DES_bs_z$nrejs[1:10]

# OrderShapeEM
ordershapeem_DES_x1_rejs <- sapply(q_levels, function(q) {
  sum(ordershape_DES_x1$fdr <= q)
})

ordershapeem_DES_x2_rejs <- sapply(q_levels, function(q) {
  sum(ordershape_DES_x2$fdr <= q)
})

# FDRreg
rdrreg_DES_rejs <- sapply(q_levels, function(q) {
  length(which(fdr_DES_theo$FDR <= q))
})

# IHW
ihw_DES_x1_rejs <- c(rejections(ihw_DES_x1_01), rejections(ihw_DES_x1_02),
                     rejections(ihw_DES_x1_03), rejections(ihw_DES_x1_04),
                     rejections(ihw_DES_x1_05), rejections(ihw_DES_x1_06),
                     rejections(ihw_DES_x1_07), rejections(ihw_DES_x1_08),
                     rejections(ihw_DES_x1_09), rejections(ihw_DES_x1_1))

ihw_DES_x2_rejs <- c(rejections(ihw_DES_x2_01), rejections(ihw_DES_x2_02),
                     rejections(ihw_DES_x2_03), rejections(ihw_DES_x2_04),
                     rejections(ihw_DES_x2_05), rejections(ihw_DES_x2_06),
                     rejections(ihw_DES_x2_07), rejections(ihw_DES_x2_08),
                     rejections(ihw_DES_x2_09), rejections(ihw_DES_x2_1))


# BH
bh_rejs <- sapply(q_levels, function(q) {
  length(which(DES_p_adjusted <= q))
})


# combine into data frame
gopa_sum_16s_DES <- data.frame(
  q = q_levels,
  chai = chai_DES_rejs,
  adapt_glm = adapt_glm_DES_ns_rej,
  # adapt_glm_bs = adapt_glm_DES_bs_rej,
  adapt_gmm_p = adapt_gmm_DES_p_ns_rej,
  # adapt_gmm_p_bs = adapt_gmm_DES_p_bs_rej,
  adapt_gmm_z = adapt_gmm_DES_z_ns_rej,
  # adapt_gmm_z_bs = adapt_gmm_DES_z_bs_rej,
  OrderShapeEM_x1 = ordershapeem_DES_x1_rejs,
  OrderShapeEM_x2 = ordershapeem_DES_x2_rejs,
  FDRreg = rdrreg_DES_rejs,
  IHW_x1 = ihw_DES_x1_rejs,
  IHW_x2 = ihw_DES_x2_rejs,
  DESeq2 = DESeq,
  BH = bh_rejs
)



###########################################################
# Plots
library(ggplot2)
library(tidyr)
library(ggrepel)
library(viridisLite)

# into long format
gopa_16s_DES_long <- gopa_sum_16s_DES %>%
  pivot_longer(
    cols = -q,
    names_to = "Method",
    values_to = "Rejections"
  ) %>%
  mutate(Method = as.character(Method))

p_gopa_16s_pcoa_DES <- plot_rejections_vs_q(gopa_16s_DES_long,
                                        title = "Number of Discoveries vs q, on DESeq2 statistics",
                                        x_breaks = sort(unique(gopa_16s_DES_long$q)))

pdf("./plots/Gopalakrishnan/16s_DES_pcoa_family_updated.pdf", width = 8, height = 4)
plot(p_gopa_16s_pcoa_DES)
dev.off()


# save(gopa_sum_16s_family, gopa_sum_16s_DES, file = "./chai_env/Gopalakrishnan_16s_tables.RData", compress = "xz")




# # into long format
# eraw_long_16s_DES <- eraw_sum_16s_DES %>%
#   pivot_longer(
#     cols = -q,
#     names_to = "Method",
#     values_to = "Rejections"
#   )
#
# eraw_chai_DES <- filter(eraw_long_16s_DES, Method == "chai")
# eraw_others_DES <- filter(eraw_long_16s_DES, Method != "chai")
#
#
# # install.packages("ggrepel")  # once
# library(ggrepel)
# library(viridisLite)
#
# # color map: viridis for others, red for chai
# methods <- unique(eraw_long_16s_DES$Method)
# others  <- setdiff(methods, "chai")
# pal     <- viridis(length(others), option = "plasma")
# method_cols <- c(setNames(pal, others), "chai" = "red")
#
# # last point (max q) for each method
# endpoints <- eraw_long_16s_DES |>
#   dplyr::arrange(Method, q) |>
#   dplyr::group_by(Method) |>
#   dplyr::slice_tail(n = 1) |>
#   dplyr::ungroup()
#
# # plot with end-of-line labels
# p <- ggplot() +
#   geom_line(data = eraw_others_DES, aes(q, Rejections, color = Method), size = 0.9) +
#   geom_point(data = eraw_others_DES, aes(q, Rejections, color = Method), size = 1.8) +
#   geom_line(data = eraw_chai_DES,   aes(q, Rejections, color = Method), size = 1.5) +
#   geom_point(data = eraw_chai_DES,  aes(q, Rejections, color = Method), size = 2.5) +
#   geom_text_repel(
#     data = endpoints,
#     aes(q, Rejections, label = Method, color = Method),
#     nudge_x = diff(range(eraw_long_16s_DES$q)) * 0.03,  # push labels a bit to the right
#     hjust = 0, direction = "y",
#     segment.size = 0.2, box.padding = 0.1, point.padding = 0.1,
#     show.legend = FALSE
#   ) +
#   scale_color_manual(values = method_cols) +
#   labs(title = "Number of Discoveries vs q", x = "q level", y = "Number of Rejections") +
#   coord_cartesian(clip = "off") +
#   expand_limits(x = max(eraw_long_16s_DES$q) + 0.01) +         # room for labels
#   theme_minimal() +
#   theme(
#     plot.margin = margin(5.5, 60, 5.5, 5.5),           # extra right margin
#     legend.position = "none"                            # labels replace legend; set to "bottom" if you want both
#   )








######################################################################################
############################ Log of normalized count #################################
######################################################################################
# Side info using log-transformed abundance table
log_count <- log1p(colSums(genera_abun))


####################################################################
############################ PLOTS #################################
####################################################################
# Wilcoxon PCoA
# Explore the selections
# Visualization related to chai results
chai_sel_pcoa <- lFDRselect(chai_family_pcoa$lFDR, 0.05, 1)

sel_otus_pcoa <- names(z_value)[chai_sel_pcoa]   # selected OTU IDs
intersect(sel_otus, sel_otus_pcoa) # Picked by both pcoa and LLM

# Tree
tr <- phy_tree(family)                 # phylo
sel_otus_pcoa <- intersect(sel_otus_pcoa, tr$tip.label)

# prune tree to selected OTUs
tr_sel_pcoa <- keep.tip(tr, sel_otus_pcoa)

# Brownian-motion phylogenetic covariance, then correlation
V_pcoa <- vcv.phylo(tr_sel_pcoa)     # covariance implied by tree branch lengths
R_phylo_pcoa <- cov2cor(V_pcoa)      # correlation matrix

# Rename OTU to taxa
tax <- as.data.frame(tax_table(family))
tax$OTU <- rownames(tax)
lab_df_pcoa <- tax %>%
  filter(OTU %in% sel_otus_pcoa) %>%
  mutate(
    raw = .data[["family"]],
    clean = clean_taxon(raw),
    label = ifelse(is.na(clean) | clean == "", OTU, paste0(clean, " (", OTU, ")"))
  )

lbl_pcoa <- setNames(make.unique(lab_df_pcoa$label), lab_df_pcoa$OTU)

rownames(R_phylo_pcoa) <- lbl_pcoa[rownames(R_phylo_pcoa)]
colnames(R_phylo_pcoa) <- lbl_pcoa[colnames(R_phylo_pcoa)]
R_phylo_pcoa[1:5, 1:5]

library(pheatmap)
pheatmap(
  R_phylo_pcoa,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "Phylogenetic correlation among selected taxa (PCoA)"
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
tip_cols <- ifelse(tr_all$tip.label %in% sel_otus_pcoa, "red", "black")

# plot with original OTU labels
plot(tr_all, show.tip.label = FALSE, cex = 0.7, tip.color = tip_cols,
     main = "Phylogenetic tree with selected families by chai using PCoA")

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

anno_pcao <- tax_all %>%
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
      OTU %in% sel_otus_pcoa & OTU %in% bh_otus ~ "both",
      OTU %in% sel_otus_pcoa & !(OTU %in% bh_otus) ~ "chai_only",
      !(OTU %in% sel_otus_pcoa) & OTU %in% bh_otus ~ "BH_only",
      TRUE ~ "neither"
    )  ) %>%
  dplyr::select(OTU, tip_label, phylum_clean, tip_color) %>%
  distinct(OTU, .keep_all = TRUE)     # 确保顺序与树 tip 一致（ggtree 依赖 label 对齐）

tr <- ape::as.phylo(tr_all)
p <- ggtree(tr, layout = "circular", size = 0.5)

p1 <- p %<+% anno_pcao

# 外圈 phylum 色块（每个 tip 一个小方块）
p2 <- p1 +
  new_scale_fill() +
  ggtreeExtra::geom_fruit(
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
  labs(title = "Phylogenetic tree with selected families by chai on wilcoxon z with PCoA") +
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

# pdf("./plots/Gopalakrishnan/16s_wilcoxz_PCoA_phylo_tree.pdf", width = 8, height = 6.5)
# plot(p4)
# dev.off()



# DESeq2 PCoA
sel_otus_DES <- rownames(res_05)[chai_DES_sel]   # selected OTU IDs
intersect(sel_otus, sel_otus_DES) # Picked by both LLM and DESeq2
# [1] "OTU_7"   "OTU_175" "OTU_95"  "OTU_646" "OTU_375" "OTU_479"

# Tree
tr <- phy_tree(family)                 # phylo
sel_otus_DES <- intersect(sel_otus_DES, tr$tip.label)

# prune tree to selected OTUs
tr_sel_DES <- keep.tip(tr, sel_otus_DES)

# Brownian-motion phylogenetic covariance, then correlation
V_DES <- vcv.phylo(tr_sel_DES)     # covariance implied by tree branch lengths
R_phylo_DES <- cov2cor(V_DES)      # correlation matrix

# Rename OTU to taxa
tax <- as.data.frame(tax_table(family))
tax$OTU <- rownames(tax)
lab_df_DES <- tax %>%
  filter(OTU %in% sel_otus_DES) %>%
  mutate(
    raw = .data[["family"]],
    clean = clean_taxon(raw),
    label = ifelse(is.na(clean) | clean == "", OTU, paste0(clean, " (", OTU, ")"))
  )

lbl_DES <- setNames(make.unique(lab_df_DES$label), lab_df_DES$OTU)

rownames(R_phylo_DES) <- lbl_DES[rownames(R_phylo_DES)]
colnames(R_phylo_DES) <- lbl_DES[colnames(R_phylo_DES)]
R_phylo_DES[1:5, 1:5]

library(pheatmap)
pheatmap(
  R_phylo_DES,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "Phylogenetic correlation among selected taxa (DESeq2 + PCoA)"
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
sel_otus_BH <- rownames(res_05)[DES_bh]   # BH selected OTU IDs
in_bh  <- tr_all$tip.label %in% sel_otus_BH

in_sel <- tr_all$tip.label %in% sel_otus_DES

tip_cols <- ifelse(in_sel & in_bh, "green",
                   ifelse(in_sel & !in_bh, "red",
                          ifelse(!in_sel & in_bh, "blue", "black")))

# plot with original OTU labels
plot(tr_all, show.tip.label = FALSE, cex = 0.7, tip.color = tip_cols,
     main = "Phylogenetic tree with selected families by chai using PCoA")

# add cleaned labels on top (so you see taxa names)
tiplabels(text = unname(lbl_all[tr_all$tip.label]) , tip = 1:length(tr_all$tip.label),
          frame = "none", adj = 0, cex = 0.7, col = tip_cols)

legend("bottomright",
       legend = c("Only in chai", "Only in DESeq2", "In both", "Neither"),
       col    = c("red", "blue", "green", "black"),
       pch    = 16,
       bty    = "n",
       cex    = 0.8,        # smaller text
       pt.cex = 0.9,        # point size
       x.intersp = 0.5,     # tighten gap between point and text
       y.intersp = 0.2,     # tighten vertical spacing between rows
       inset = 0.02)


# Circular tree
library(dplyr)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

rank_to_show <- "family"

anno_DES <- tax_all %>%
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
      OTU %in% sel_otus_DES & OTU %in% sel_otus_BH ~ "both",
      OTU %in% sel_otus_DES & !(OTU %in% sel_otus_BH) ~ "chai_only",
      !(OTU %in% sel_otus_DES) & OTU %in% sel_otus_BH ~ "BH_only",
      TRUE ~ "neither"
    )  ) %>%
  dplyr::select(OTU, tip_label, phylum_clean, tip_color) %>%
  distinct(OTU, .keep_all = TRUE)     # 确保顺序与树 tip 一致（ggtree 依赖 label 对齐）

tr <- ape::as.phylo(tr_all)
p <- ggtree(tr, layout = "circular", size = 0.5)

p1 <- p %<+% anno_DES

# 外圈 phylum 色块（每个 tip 一个小方块）
p2 <- p1 +
  new_scale_fill() +
  ggtreeExtra::geom_fruit(
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
    offset = 0.02,          # 文字离色块距离：想更远就调大
    show.legend = FALSE
  ) +
  scale_color_manual(
    name   = "Tip label group",
    values = c(chai_only="red", BH_only="blue", both="green", neither="black"),
    breaks = c("chai_only", "BH_only", "both", "neither"),
    labels = c("Only in chai", "Only in DESeq2", "In both", "Neither"),
    drop   = FALSE
  ) +
  labs(title = "Phylogenetic tree with selected families by chai on DESeq2 z with PCoA") +
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

# pdf("./plots/Gopalakrishnan/16s_DESeq2_PCoA_phylo_tree.pdf", width = 8, height = 6.5)
# plot(p4)
# dev.off()
