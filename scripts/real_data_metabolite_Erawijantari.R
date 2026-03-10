# Erawijantari's data - Gastric cancer vs. control - 2020
library(dplyr)
library(vegan)
library(ggplot2)
library(ggpubr)
library(adaptMT)
library(AdaPTGMM)
library(splines)
library(splines2)
library(FDRreg)
library(coin)
library(IHW)

setwd("C:/Users/zwang26/OneDrive - UTHealth Houston/conditionalGaussian")
# source("chai.R")
source("color_helper.R")

# load("./chai_env/eraw_tables.RData")

########################## Data from Github #################################

mtb <- read.delim("./Erawijantari/mtb.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
# 96 x 524
genera_counts <- read.delim("./Erawijantari/genera.counts.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
# 96 x 10527
metadata <- read.delim("./Erawijantari/metadata.tsv", header = TRUE, stringsAsFactors = FALSE)


table(metadata$Study.Group)  # With metabolite profiles
# Gastrectomy 42 vs Healthy 54

table(rownames(mtb) == rownames(genera_counts)) # The subjects in mtb is matching the subjects in genera_counts

####################### genera_counts ############################
##################### Filter out < 10% ###########################
# Filter out the 'low abundance' across all subjects
# Keep genera that are non-zero in at least 10% of subjects
present_counts <- colSums(genera_counts > 0)
genera_counts_filtered <- genera_counts[, present_counts >= (nrow(genera_counts) * 0.1)]
# 96 x 5163

# Get the abundance table - genera_abun
row_sums <- rowSums(genera_counts_filtered)
genera_abun <- sweep(genera_counts_filtered, 1, row_sums, FUN="/")
# 96 x 5163

# Re-order the metadata to match the order of genera.
metadata_reorder <- metadata[match(rownames(genera_abun), metadata$Subject), ]


######################################################################################
######################## Correlation of metabolite data ##############################
######################################################################################

######################## Metabolite ##############################
##################### Filter out < 50% ###########################
present_counts <- colSums(mtb > 0)
mtb_filtered <- mtb[, present_counts >= (nrow(mtb) * 0.5)]
dim(mtb_filtered)
# 96 x 144

# Get top 50 importance of metabolites
mtb_norm <- mtb_filtered / rowSums(mtb_filtered)
metabolite_means <- colMeans(mtb_norm, na.rm = TRUE)
top50_metabolites <- names(sort(metabolite_means, decreasing = TRUE))[1:50]
mtb_top50 <- mtb_norm[, top50_metabolites]


# Check the order
# making sure that rownames of meta and of abun are matching
table(rownames(genera_abun) == rownames(mtb_top50))
table(rownames(genera_abun) == metadata_reorder$Subject)


################# Canonical Correlation #####################
cancor <- numeric(ncol(genera_abun))
ycoef <- matrix(0, nrow = ncol(genera_abun), ncol = ncol(mtb_top50))

for (i in 1:ncol(genera_abun)) {
  temp <- cancor(genera_abun[, i], mtb_top50)
  cancor[i]  <- temp$cor
  ycoef[i, ] <- temp$ycoef[, 1]
}

# Get the outcome - outcome
outcome <- as.factor(metadata_reorder$Study.Group)


##################### p-value/z-value ###########################
p_z_values <- sapply(colnames(genera_abun), function(feature) {
  features <- genera_abun[, feature, drop = FALSE]
  wilcox <- wilcox_test(features[,1] ~ outcome, distribution = "asymptotic")

  p_value <- as.numeric(coin::pvalue(wilcox))
  z_value <- statistic(wilcox)

  return(c(p_value = p_value, z_value = z_value))
}, simplify = "matrix")
p_z_values <- t(p_z_values)

p_value <- p_z_values[,1]
z_value <- p_z_values[,2]

# p_cal <- 2 * (1 - pnorm(abs(z_value)))

##################### BH ###########################
p_adjusted <- p.adjust(p_value, method = "BH")
bh <- which(p_adjusted <= 0.05)
length(bh) #98


################## Visualization ########################
df <- data.frame(z = z_value, x = cancor)
# z
ggplot(df, aes(x=z)) +
  geom_histogram(aes(y=..density..), binwidth=0.25, color="black", fill="lightblue") +
  geom_density() +
  stat_function(fun=dnorm, args=list(mean=0, sd=1), linetype="dashed") +
  labs(title="Histogram of z-value with N(0,1)", x="z", y="Density") +
  theme_minimal()

# x - Canonical Correlation
ggplot(df, aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 40, color = "black", fill = "lightblue") +
  geom_density() +
  labs(
    title = paste0("Histogram of X - Canonical Correlation"),
    x = "X", y = "Density"
  ) +
  theme_minimal()

# z vs. x
ggplot(df, aes(x = z, y = x)) +
  geom_point(alpha = 0.6) +
  labs(x = "z-value", y = "Canonical Correlation") +
  theme_minimal()

# p
ggplot(data.frame(p_value), aes(p_value)) +
  geom_histogram(bins=40, color="black", fill="grey80") +
  labs(title="Histogram of two-sided p-values", x="p", y="Count") +
  theme_minimal()

####################################################################
######################## Explore data ##############################
####################################################################
# alpha diversity
alpha_df <- tibble(
  sample = rownames(genera_abun),
  group = outcome,
  shannon = vegan::diversity(genera_abun, index = "shannon"),
  simpson = vegan::diversity(genera_abun, index = "simpson")
)

alpha_df$group <- factor(alpha_df$group, levels = c("Gastrectomy", "Healthy"))

cols <- c("Gastrectomy" = "#C8A100",
          "Healthy" = "#2C7FB8")

a_p <- wilcox.test(shannon ~ group, data = alpha_df)$p.value
a_lab <- sprintf("p-value = %.3g", a_p)


pa <- ggplot(alpha_df, aes(x = group, y = shannon, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.25) +
  geom_jitter(aes(color = group), width = 0.15, alpha = 0.6) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  labs(x = NULL, y = "Shannon", title = "Alpha diversity (Shannon)") +
  annotate("text",
           x = 1.5, y = max(alpha_df$shannon, na.rm = TRUE)*1.2,
           label = a_lab,
           size = 4, color = "black") +
  theme_bw() +
  theme(legend.position = "none")


# beta diversity
bc <- vegan::vegdist(genera_abun, method = "bray")

# PCoA
pcoa_res <- ape::pcoa(bc)
coords <- as.data.frame(pcoa_res$vectors[, 1:2, drop = FALSE])
colnames(coords) <- c("PCoA1", "PCoA2")
coords$sample <- rownames(coords)
coords$group <- outcome

# variance explained
rel_eig <- pcoa_res$values$Relative_eig
ve1 <- round(100*rel_eig[1], 1)
ve2 <- round(100*rel_eig[2], 1)

# PERMANVOA
ad <- adonis2(bc ~ outcome, permutations = 9999)

ad_r2 <- ad$R2[1]
ad_p <- ad$`Pr(>F)`[1]
ad_lab <- sprintf("PERMANOVA: R2 = %.3f, p = %.3g", ad_r2, ad_p)

coords$group <- factor(coords$group, levels = c("Gastrectomy", "Healthy"))

ppcoa <- ggplot(coords, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = group), size = 2.5, alpha = 0.85) +
  stat_ellipse(aes(color = group, fill = group),
               level = 0.95, type = "t", linewidth = 1) +
  annotate("text",
           x = min(coords$PCoA1, na.rm = TRUE),
           y = max(coords$PCoA2, na.rm = TRUE),
           label = ad_lab,
           hjust = 0, vjust = 1,
           size = 4, color = "black") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(
    title = "PCoA on Bray-Curtis distance",
    x = paste0("PCoA1 (", ve1, "%)"),
    y = paste0("PCoA2 (", ve2, "%)")
  ) +
  theme_bw()

library(patchwork)
p_ab <- pa + ppcoa +
  plot_layout(ncol = 2, widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "A")

print(p_ab)

# pdf("./plots/Erawijantari/Erawijantari_a_b_diversity.pdf", width = 10, height = 4)
# plot(p_ab)
# dev.off()

# cancor vs. z
df <- data.frame(z = z_value, x = cancor)
p_cca_z <- ggplot(df, aes(x = z, y = x)) +
            geom_point(alpha = 0.6) +
            labs(x = "z-value", y = "Canonical Correlation") +
            theme_minimal()
# pdf("./plots/Erawijantari/Erawijantari_cca_z.pdf", width = 5, height = 4)
# plot(p_cca_z)
# dev.off()

############################ Fit model ############################
################## chai #######################
set.seed(123)
chai_cca <- chai(z_value, as.matrix(cancor, ncol = 1), R = 100)

# select <- lFDRselect(chai_cca, q = 0.05)
# length(lFDRselect(chai_cca, q = 0.05))

chai_sel <- lFDRselect(chai_cca$lFDR, 0.05, 1)

length(lFDRselect(chai_cca$lFDR, 0.01, 1)) # 63
length(lFDRselect(chai_cca$lFDR, 0.05, 1)) # 208
length(lFDRselect(chai_cca$lFDR, 0.1, 1))  # 379


label_direction <- function(outcome, z_value, selected_idx) {
  first_level <- levels(outcome)[1]
  lab_hi <- paste0("higher_in_", first_level)
  lab_lo <- paste0("lower_in_", first_level)
  data.frame(
    feature   = names(z_value)[selected_idx],
    z_value   = unname(z_value[selected_idx]),
    direction = ifelse(z_value[selected_idx] > 0, lab_hi, lab_lo),
    stringsAsFactors = FALSE
  )
}

chai_sel_dir <- label_direction(outcome, z_value, chai_sel)


# Visualization
df <- data.frame(z = z_value, x = cancor) %>%
  dplyr::mutate(
    grp = dplyr::case_when(
      dplyr::row_number() %in% bh & dplyr::row_number() %in% chai_sel ~ "both",
      dplyr::row_number() %in% bh ~ "BH_only",
      dplyr::row_number() %in% chai_sel ~ "chai_only",
      TRUE ~ "neither"
    )
  )

cca_sel <- ggplot(df, aes(x = z, y = x, color = grp)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(
    values = c(
    chai_only  = "red",
    BH_only  = "blue",
    both    = "green",
    neither = "black"
    ),
    breaks = c("chai_only", "BH_only", "both", "neither"),
    labels = c("Only in chai", "Only in BH", "In both", "Neither")
    ) +
  labs(title = "z-statistics vs CCA coefficient",
       x = "z-value",
       y = "Canonical Correlation",
       color = "Group") +
  theme_minimal()

# pdf("./plots/Erawijantari/Erawijantari_cca_z.pdf", width = 8, height = 6)
# plot(cca_sel)
# dev.off()

# Check the feature names
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

# feature names (208 of them)
feat <- names(z_value)[chai_sel]

get_rank <- function(x, rank_prefix) {
  # rank_prefix like "d", "p", "c", "o", "f", "g"
  # grabs text after "{rank}__" up to the next ".<letter>__" or end of string
  str_match(x, paste0(rank_prefix, "__(.*?)(?=\\.[a-z]__|$)"))[, 2]
}

chai_tax_tbl <- chai_sel_dir %>%
  as_tibble() %>%
  mutate(
    domain = get_rank(feature, "d"),
    phylum = get_rank(feature, "p"),
    class  = get_rank(feature, "c"),
    order  = get_rank(feature, "o"),
    family = get_rank(feature, "f"),
    genus  = get_rank(feature, "g")
  ) %>%
  mutate(across(
    c(domain, phylum, class, order, family, genus),
    ~ .x %>%
      str_replace_all("\\[|\\]", "") %>%
      str_squish()
  )) %>%
  # keep your original columns + taxonomy columns (order them nicely)
  dplyr::select(feature, z_value, direction, domain, phylum, class, order, family, genus)

# Collapse the *_A, *_C, and so on
chai_tax_tbl2 <- chai_tax_tbl %>%
  mutate(
    phylum_collapsed = sub("_(A|B|C|D|E|F|G|H|I|J|K)$", "", phylum)
  )

chai_tax_tbl %>% dplyr::count(phylum, sort=TRUE)
chai_tax_tbl2 %>% dplyr::count(phylum_collapsed, sort = TRUE)
chai_tax_tbl2 %>%
  dplyr::count(phylum_collapsed, direction, sort = TRUE)


########################### Likert plot ###########################
library(dplyr)
library(ggplot2)

# 1. Prepare the data
plot_data <- chai_tax_tbl2 %>%
  dplyr::count(phylum_collapsed, direction) %>%
  # Make 'lower' counts negative so they extend to the left of zero
  mutate(count_adj = ifelse(direction == "lower_in_Gastrectomy", -n, n))

# 2. Determine the order (sorted by most 'higher' counts)
phylum_order <- chai_tax_tbl2 %>%
  dplyr::count(phylum_collapsed, sort=TRUE) %>%
  pull(phylum_collapsed)

# 3. Create the plot
ggplot(plot_data, aes(y = factor(phylum_collapsed, levels = phylum_order),
                      x = count_adj,
                      fill = direction)) +
  geom_col() +
  # Use scale_x_continuous to show absolute values (removing the minus signs)
  scale_x_continuous(labels = abs) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  theme_minimal() +
  labs(
    title = "Phylum Distribution: Lower vs. Higher",
    x = "Genera count",
    y = "Phylum",
    fill = "Direction"
  ) +
  scale_fill_manual(values = c("higher_in_Gastrectomy" = "#80cdc1", "lower_in_Gastrectomy" = "#dfc27d"))

#################################################################################


# BH found:
bh_feat <- names(z_value)[bh]
tax_bh <- tibble(feature = bh_feat) %>%
  mutate(
    domain = get_rank(feature, "d"),
    phylum = get_rank(feature, "p"),
    class  = get_rank(feature, "c"),
    order  = get_rank(feature, "o"),
    family = get_rank(feature, "f"),
    genus  = get_rank(feature, "g")
  ) %>%
  mutate(across(c(domain, phylum, class, order, family, genus),
                ~ .x %>% str_replace_all("\\[|\\]", "") %>% str_squish()))

# Collapse the *_A, *_C, and so on
tax_bh2 <- tax_bh %>%
  mutate(
    phylum_collapsed = sub("_(A|B|C|D|E|F|G|H|I|J|K)$", "", phylum)
  )

tax_bh2 %>% dplyr::count(phylum, sort=TRUE)
tax_bh2 %>% dplyr::count(phylum_collapsed, sort = TRUE)


# species PCA
spe_met_pca <- prcomp(ycoef, center = TRUE, scale. = TRUE)

# Build plotting dataframe (PC scores for each species/row)
df_pca <- data.frame(
  idx = seq(nrow(ycoef)),                # original row number in ycoef
  PC1 = spe_met_pca$x[, 1],
  PC2 = spe_met_pca$x[, 2],
  stringsAsFactors = FALSE
)

in_c <- df_pca$idx %in% chai_sel
in_b <- df_pca$idx %in% bh

df_pca$group_color <- ifelse(in_c & in_b, "both",
                             ifelse(in_c & !in_b, "chai_only",
                                    ifelse(!in_c & in_b, "BH_only", "neither")))
# Plot
library(ggplot2)
ggplot(df_pca, aes(x = PC1, y = PC2, color = group_color)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_color_manual(
    values = c(
      chai_only  = "red",
      BH_only  = "blue",
      both    = "green",
      neither = "black"
    ),
    breaks = c("chai_only", "BH_only", "both", "neither"),
    labels = c("only in chai", "only in BH", "both", "neither")
  ) +
  theme_minimal() +
  labs(
    title = "PCA of CCA coefficient between species vs. metabolites",
    x = paste0("PC1 (", round(100 * summary(spe_met_pca)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(spe_met_pca)$importance[2, 2], 1), "%)"),
    color = "Group"
  )


df_show <- subset(df_pca, group_color != "neither")
ggplot(df_show, aes(x = PC1, y = PC2, color = group_color)) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_manual(
    values = c(
      chai_only = "red",
      BH_only = "blue",
      both   = "green"
    ),
    breaks = c("chai_only", "BH_only", "both"),
    labels = c("only in chai", "only in BH", "in both")
  ) +
  theme_minimal() +
  labs(
    title = "PCA of CCA coefficient between species vs. metabolites (only selected points)",
    x = paste0("PC1 (", round(100 * summary(spe_met_pca)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(spe_met_pca)$importance[2, 2], 1), "%)"),
    color = "Group"
  )

################## Adapt #######################
alphas <- seq(0.01, 0.10, by = 0.01)

# GLM
# ns
x <- data.frame(cancor)
formula <- c(paste0("ns(cancor, df = ",2:8,")"))
glm_ns <- adapt_glm(x, pvals = p_value, alphas=alphas,
                    pi_formulas = formula, mu_formulas = formula)
# alpha = 0.1: FDPhat 0.0992, Number of Rej. 359

# # ispline
# x <- data.frame(cancor)
# formulas <- c(paste0("iSpline(cancor, df = ",2:8,")"))
# glm_ispline <- adapt_glm(x, pvals = p_value,  alphas=alphas,
#                          pi_formulas = formulas, mu_formulas = formulas)
# # alpha = 0.1: FDPhat 0.0974, Number of Rej. 349
# # alpha = 0.09: FDPhat 0.0895, Number of Rej. 324
#
# # bs
# formula_bs <- c(paste0("bs(cancor, df = ",2:8,")"))
# glm_bs <- adapt_glm(x, pvals = p_value,  alphas=alphas,
#                          pi_formulas = formula_bs, mu_formulas = formula_bs)
# # alpha = 0.1: FDPhat 0.0997, Number of Rej. 371
# # alpha = 0.09: FDPhat 0.0886, Number of Rej. 316

# GMM - p
# ns
gmm_ns_p <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
                      beta_formulas = formula)
# alpha = 0.1: FDPhat 0.0983, Number of Rej. 300
# alpha = 0.05: FDPhat 0.0495, Number of Rej. 212


# # ispline
# gmm_ispline_p <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
#                            beta_formulas = formulas)
# # alpha = 0.1: FDPhat 0.0997, Number of Rej. 311
# # alpha = 0.05: FDPhat 0.0491, Number of Rej. 224
#
# # bs
# gmm_bs_p <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
#                       beta_formulas = formula_bs)
#
# # alpha = 0.1: FDPhat 0.0984, Number of Rej. 305
# # alpha = 0.05: FDPhat 0.05, Number of Rej. 230


# GMM - Z
gmm_ns_z <- adapt_gmm(x, z = z_value, alphas=alphas,
                      beta_formulas = formula, testing = "two_sided")
# alpha = 0.1: FDPhat 0.0997, Number of Rej. 386
# alpha = 0.05: FDPhat 0.0482, Number of Rej. 197

# gmm_bs_z <- adapt_gmm(x, z = z_value, alphas=alphas,
#                            beta_formulas = formula_bs, testing = "two_sided")
# # alpha = 0.1: FDPhat 0.0992, Number of Rej. 388
# # alpha = 0.05: FDPhat 0.0477, Number of Rej. 199

################## OrderShapeEM #######################

require(OrderShapeEM)
orderfdr <- OrderShapeEM(pvals = p_value, order.var = cancor,
                         OrderShapeEM.control(trace = TRUE))
sum(orderfdr$fdr <= 0.05)
# 130

######################## FDRreg ############################
set.seed(123)
fdr_theo <- FDRreg(z_value, as.matrix(cancor), nulltype = 'theoretical')
length(which(fdr_theo$FDR <= 0.05))
# 160

fdr_emp <- FDRreg(z_value, as.matrix(cancor), nulltype = 'empirical')
length(which(fdr_emp$FDR <= 0.05))
# 8

######################## IHW ############################
ihw_01 <- ihw(pvalues = p_value, covariates = cancor, alpha = 0.01)
ihw_02 <- ihw(pvalues = p_value, covariates = cancor, alpha = 0.02)
ihw_03 <- ihw(pvalues = p_value, covariates = cancor, alpha = 0.03)
ihw_04 <- ihw(pvalues = p_value, covariates = cancor, alpha = 0.04)
ihw_05 <- ihw(pvalues = p_value, covariates = cancor, alpha = 0.05)
ihw_06 <- ihw(pvalues = p_value, covariates = cancor, alpha = 0.06)
ihw_07 <- ihw(pvalues = p_value, covariates = cancor, alpha = 0.07)
ihw_08 <- ihw(pvalues = p_value, covariates = cancor, alpha = 0.08)
ihw_09 <- ihw(pvalues = p_value, covariates = cancor, alpha = 0.09)
ihw_1 <- ihw(pvalues = p_value, covariates = cancor, alpha = 0.1)

rejections(ihw_05)  # 123
rejections(ihw_06)  # 146
rejections(ihw_07)  # 167
rejections(ihw_08)  # 180
rejections(ihw_09)  # 205
rejections(ihw_1)   # 241

###########################################################
# Summarize into a different q level table
q_levels <- alphas

# chai
chai_rejs <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_cca$lFDR, q, 1))
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
ordershapeem_rejs <- sapply(q_levels, function(q) {
  sum(orderfdr$fdr <= q)
})

# FDRreg
rdrreg_rejs <- sapply(q_levels, function(q) {
  length(which(fdr_theo$FDR <= q))
})

# IHW
ihw_rejs <- c(rejections(ihw_01), rejections(ihw_02), rejections(ihw_03),
              rejections(ihw_04), rejections(ihw_05), rejections(ihw_06),
              rejections(ihw_07), rejections(ihw_08), rejections(ihw_09),
              rejections(ihw_1))
# [1]  36  61  72 110 123 146 167 180 205 241

# BH
bh_rejs <- sapply(q_levels, function(q) {
  length(which(p_adjusted <= q))
})


# combine into data frame
eraw_sum <- data.frame(
  q = q_levels,
  chai = chai_rejs,
  adapt_glm = adapt_glm_ns_rej,
  # adapt_glm_bs = adapt_glm_bs_rej,
  adapt_gmm_p = adapt_gmm_p_ns_rej,
  # adapt_gmm_p_bs = adapt_gmm_p_bs_rej,
  adapt_gmm_z = adapt_gmm_z_ns_rej,
  # adapt_gmm_z_bs = adapt_gmm_z_bs_rej,
  OrderShapeEM = ordershapeem_rejs,
  FDRreg = rdrreg_rejs,
  IHW = ihw_rejs,
  BH = bh_rejs
)



###########################################################
# Plots
# Consistant color

library(ggplot2)
library(tidyr)
library(ggrepel)
library(viridisLite)

# into long format
eraw_long <- eraw_sum %>%
  pivot_longer(
    cols = -q,
    names_to = "Method",
    values_to = "Rejections"
  ) %>%
  mutate(Method = as.character(Method))
#
# order <- c("chai", "adapt_glm", "adapt_gmm_z", "adapt_gmm_p",
#            "FDRreg", "OrderShapeEM", "IHW", "BH")
#
# eraw_long$Method <- factor(eraw_long$Method, levels = order)

plot_rejections_vs_q <- function(eraw_long, title = "Number of Discoveries vs q",
                                 K_variants = 12, x_breaks = NULL) {

  ord <- method_order(eraw_long$Method)
  eraw_long <- eraw_long %>% mutate(Method = factor(Method, levels = ord))
  pal <- make_method_palette(levels(eraw_long$Method), K = K_variants)

  endpoints <- eraw_long %>%
    group_by(Method) %>%
    filter(q == max(q, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    mutate(label = as.character(Method))

  x_rng   <- range(eraw_long$q, na.rm = TRUE)
  x_nudge <- diff(x_rng) * 0.03

  p <- ggplot(eraw_long, aes(x = q, y = Rejections, group = Method, color = Method)) +
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

p_eraw_cca <- plot_rejections_vs_q(eraw_long,
                                   title = "Number of Discoveries vs q, with metabolites CCA score",
                                   x_breaks = sort(unique(eraw_long$q)))


#
# eraw_chai <- filter(eraw_long, Method == "chai")
# eraw_others <- filter(eraw_long, Method != "chai")
#
# # color map: viridis for others, red for chai
# methods <- unique(eraw_long$Method)
# others  <- setdiff(methods, "chai")
# pal     <- viridis(length(others), option = "plasma")
# method_cols <- c(setNames(pal, others), "chai" = "red")
#
# # last point (max q) for each method
# endpoints <- eraw_long |>
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
#     nudge_x = diff(range(eraw_long$q)) * 0.03,  # push labels a bit to the right
#     hjust = 0, direction = "y",
#     segment.size = 0.2, box.padding = 0.1, point.padding = 0.1,
#     show.legend = FALSE
#   ) +
#   scale_color_manual(values = method_cols) +
#   labs(title = "Number of Discoveries vs q", x = "q level", y = "Number of Rejections") +
#   coord_cartesian(clip = "off") +
#   expand_limits(x = max(eraw_long$q) + 0.01) +         # room for labels
#   theme_minimal() +
#   theme(
#     plot.margin = margin(5.5, 60, 5.5, 5.5),           # extra right margin
#     legend.position = "none"                            # labels replace legend; set to "bottom" if you want both
#   )


# ggsave("C:/Users/zwang26/OneDrive - UTHealth Houston/conditionalGaussian/plots/Simulations/SimZ_fdr_power_updated.png", plot = combo, width = 8, height = 12, dpi = 300)


pdf("./plots/Erawijantari/Erawijantari_cca_updated.pdf", width = 8, height = 4)
plot(p_eraw_cca)
dev.off()



######################################################################################
######################## Correlation of metabolite data ##############################
######################################################################################

# ######################## KEY Metabolite from paper ##############################
#
# key <- c("C04483", "C00695", "C00245", "C01921")
# key_metabolite <- grepl(paste(key, collapse = "|"), colnames(mtb))
# key_mtb <- mtb[, key_metabolite, drop = FALSE]
#
# ################# Canonical Correlation #####################
# cancor_key <- numeric(ncol(genera_abun))
# for (i in 1:ncol(genera_abun)) {
#   temp <- cancor(genera_abun[, i], key_mtb)
#   cancor_key[i] <- temp$cor
# }
#
#
# ################## Visualization ########################
# df <- data.frame(z = z_value, x = cancor_key)
# # z
# ggplot(df, aes(x=z)) +
#   geom_histogram(aes(y=..density..), binwidth=0.25, color="black", fill="lightblue") +
#   geom_density() +
#   stat_function(fun=dnorm, args=list(mean=0, sd=1), linetype="dashed") +
#   labs(title="Histogram of z-value with N(0,1)", x="z", y="Density") +
#   theme_minimal()
#
# # x - Canonical Correlation
# ggplot(df, aes(x = x)) +
#   geom_histogram(aes(y = ..density..), bins = 40, color = "black", fill = "lightblue") +
#   geom_density() +
#   labs(
#     title = paste0("Histogram of X - Canonical Correlation"),
#     x = "X", y = "Density"
#   ) +
#   theme_minimal()
#
# # z vs. x
# ggplot(df, aes(x = z, y = x)) +
#   geom_point(alpha = 0.6) +
#   labs(x = "z-value", y = "Canonical Correlation") +
#   theme_minimal()
#
# # p
# ggplot(data.frame(p_value), aes(p_value)) +
#   geom_histogram(bins=40, color="black", fill="grey80") +
#   labs(title="Histogram of two-sided p-values", x="p", y="Count") +
#   theme_minimal()
#
# ############################ Fit model ############################
# ################## chai #######################
# chai_cca_keymbt <- chai(z_value, as.matrix(cancor_key, ncol = 1), M = 100)
# length(lFDRselect(chai_cca_keymbt$lFDR, 0.05, 1)) # 163
# length(lFDRselect(chai_cca_keymbt$lFDR, 0.1, 1))  # 317
#
#
# ################## Adapt #######################
# alphas <- seq(0.01,0.1, by = 0.01)
#
# # GLM
# # ns
# x <- data.frame(cancor_key)
# formula <- c(paste0("ns(cancor_key, df = ",2:8,")"))
# glm_ns_keymbt <- adapt_glm(x, pvals = p_value, alphas=alphas,
#                     pi_formulas = formula, mu_formulas = formula)
# # alpha = 0.1: FDPhat 0.0979, Number of Rej. 286
#
# # ispline
# x <- data.frame(cancor_key)
# formulas <- c(paste0("iSpline(cancor_key, df = ",2:8,")"))
# glm_ispline_keymbt <- adapt_glm(x, pvals = p_value,  alphas=alphas,
#                          pi_formulas = formulas, mu_formulas = formulas)
#
#
# # bs
# formula_bs <- c(paste0("bs(cancor_key, df = ",2:8,")"))
# glm_bs_keymbt <- adapt_glm(x, pvals = p_value,  alphas=alphas,
#                     pi_formulas = formula_bs, mu_formulas = formula_bs)
# # alpha = 0.1: FDPhat 0.0996, Number of Rej. 281
#
# # GMM - p
# # ns
# gmm_ns_p_keymbt <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
#                       beta_formulas = formula)
# # alpha = 0.1: FDPhat 0.0994, Number of Rej. 357
# # alpha = 0.05: FDPhat 0.0489, Number of Rej. 225
#
#
# # ispline
# gmm_ispline_p_keymbt <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
#                            beta_formulas = formulas)
#
# # bs
# gmm_bs_p_keymbt <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
#                       beta_formulas = formula_bs)
# # alpha = 0.1: FDPhat 0.0994, Number of Rej. 362
# # alpha = 0.05: FDPhat 0.0482, Number of Rej. 228
#
#
#
# # GMM - Z
# gmm_ns_z_keymbt <- adapt_gmm(x, z = z_value, alphas=alphas,
#                       beta_formulas = formula, testing = "two_sided")
# # alpha = 0.1: FDPhat 0.0988, Number of Rej. 425
# # alpha = 0.05: FDPhat 0.0495, Number of Rej. 283
#
# gmm_bs_z_keymbt <- adapt_gmm(x, z = z_value, alphas=alphas,
#                       beta_formulas = formula_bs, testing = "two_sided")
# # alpha = 0.1: FDPhat 0.0997, Number of Rej. 396
# # alpha = 0.05: FDPhat 0.0493, Number of Rej. 274
#
# ################## OrderShapeEM #######################
#
# require(OrderShapeEM)
# orderfdr_keymbt <- OrderShapeEM(pvals = p_value, order.var = cancor_key,
#                          OrderShapeEM.control(trace = TRUE))
# sum(orderfdr_keymbt$fdr <= 0.05)
# # 133
#
# ######################## FDRreg ############################
# fdr_theo_keymbt <- FDRreg(z_value, cancor_key, nulltype = 'theoretical')
# length(which(fdr_theo_keymbt$FDR <= 0.05))
# # 312
#
# fdr_emp_keymbt <- FDRreg(z_value, cancor_key, nulltype = 'empirical')
# length(which(fdr_emp_keymbt$FDR <= 0.05))
# # 8
#
# ###########################################################
# # Summarize into a different q level table
# q_levels <- seq(0.05, 0.10, by = 0.01)
#
# # chai
# chai_keymbt_rejs <- sapply(q_levels, function(q) {
#   length(lFDRselect(chai_cca_keymbt$lFDR, q, 1))
# })
#
# # adapt_glm
# adapt_glm_ns_keymbt_rej <- glm_ns_keymbt$nrejs[5:10]
# adapt_glm_bs_keymbt_rej <- glm_bs_keymbt$nrejs[5:10]
#
# # adapt_gmm_p
# adapt_gmm_p_ns_keymbt_rej <- gmm_ns_p_keymbt$nrejs[5:10]
# adapt_gmm_p_bs_keymbt_rej <- gmm_bs_p_keymbt$nrejs[5:10]
#
# # adapt_gmm_z
# adapt_gmm_z_ns_keymbt_rej <- gmm_ns_z_keymbt$nrejs[5:10]
# adapt_gmm_z_bs_keymbt_rej <- gmm_bs_z_keymbt$nrejs[5:10]
#
# # OrderShapeEM
# ordershapeem_keymbt_rejs <- sapply(q_levels, function(q) {
#   sum(orderfdr_keymbt$fdr <= q)
# })
#
# # FDRreg
# rdrreg_keymbt_rejs <- sapply(q_levels, function(q) {
#   length(which(fdr_theo_keymbt$FDR <= q))
# })
#
#
# # combine into data frame
# eraw_sum_keymbt <- data.frame(
#   q = q_levels,
#   chai = chai_keymbt_rejs,
#   adapt_glm_ns = adapt_glm_ns_keymbt_rej,
#   adapt_glm_bs = adapt_glm_bs_keymbt_rej,
#   adapt_gmm_p_ns = adapt_gmm_p_ns_keymbt_rej,
#   adapt_gmm_p_bs = adapt_gmm_p_bs_keymbt_rej,
#   adapt_gmm_z_ns = adapt_gmm_z_ns_keymbt_rej,
#   adapt_gmm_z_bs = adapt_gmm_z_bs_keymbt_rej,
#   ordershapeem = ordershapeem_keymbt_rejs,
#   fdrreg = rdrreg_keymbt_rejs
# )
#
#
# ###########################################################
# # Plots
# library(ggplot2)
# library(tidyr)
#
# # into long format
# eraw_long_keymbt <- eraw_sum_keymbt %>%
#   pivot_longer(
#     cols = -q,
#     names_to = "Method",
#     values_to = "Rejections"
#   )
#
# eraw_chai_keymbt <- filter(eraw_long_keymbt, Method == "chai")
# eraw_others_keymbt <- filter(eraw_long_keymbt, Method != "chai")
#
# # plot
# p_keymbt <- ggplot() +
#   # other methods with gradient palette
#   geom_line(data = eraw_others_keymbt, aes(x = q, y = Rejections, color = Method), size = 0.9) +
#   geom_point(data = eraw_others_keymbt, aes(x = q, y = Rejections, color = Method), size = 1.8) +
#   scale_color_viridis_d(option = "plasma") +
#
#   # chai method highlighted in red
#   geom_line(data = eraw_chai_keymbt, aes(x = q, y = Rejections), color = "red", size = 1.5) +
#   geom_point(data = eraw_chai_keymbt, aes(x = q, y = Rejections), color = "red", size = 2.5) +
#
#   labs(
#     title = "Number of Discoveries vs q",
#     x = "q level",
#     y = "Number of Rejections"
#   ) +
#   theme_minimal() +
#   theme(legend.position = "bottom")



######################################################################################
############################ Log of normalized count #################################
######################################################################################
# Side info using log-transformed abundance table
log_count <- log1p(colSums(genera_abun))

################## Visualization ########################
df <- data.frame(z = z_value, x = log_count)
# z
ggplot(df, aes(x=z)) +
  geom_histogram(aes(y=..density..), binwidth=0.25, color="black", fill="lightblue") +
  geom_density() +
  stat_function(fun=dnorm, args=list(mean=0, sd=1), linetype="dashed") +
  labs(title="Histogram of z-value with N(0,1)", x="z", y="Density") +
  theme_minimal()

# x - log of normalized counts
ggplot(df, aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 40, color = "black", fill = "lightblue") +
  geom_density() +
  labs(
    title = paste0("Histogram of X - log of normalized counts"),
    x = "X", y = "Density"
  ) +
  theme_minimal()

# z vs. x
ggplot(df, aes(x = z, y = x)) +
  geom_point(alpha = 0.6) +
  labs(x = "z-value", y = "log of normalized counts") +
  theme_minimal()

# p
ggplot(data.frame(p_value), aes(p_value)) +
  geom_histogram(bins=40, color="black", fill="grey80") +
  labs(title="Histogram of two-sided p-values", x="p", y="Count") +
  theme_minimal()

############################ Fit model ############################
################## chai #######################
set.seed(123)
chai_logcount <- chai(z_value, as.matrix(log_count, ncol = 1), M = 100)
# lFDRselect(chai_logcount, q = 0.05)
length(lFDRselect(chai_logcount$lFDR, 0.05, 1)) # 249
length(lFDRselect(chai_logcount$lFDR, 0.1, 1))  # 404


################## Adapt #######################
alphas <- seq(0.01, 0.10, by = 0.01)

# GLM
# ns
x <- data.frame(log_count)
formula <- c(paste0("ns(log_count, df = ",2:8,")"))
glm_ns_logcount <- adapt_glm(x, pvals = p_value, alphas=alphas,
                    pi_formulas = formula, mu_formulas = formula)
# alpha = 0.1: FDPhat 0.0992, Number of Rej. 373
# alpha = 0.09: FDPhat 0.089, Number of Rej. 292
# alpha = 0.08: FDPhat 0.0773, Number of Rej. 220

# # ispline
# formula_ispline <- c(paste0("iSpline(log_count, df = ",2:8,")"))
# glm_ispline_logcount <- adapt_glm(x, pvals = p_value,  alphas=alphas,
#                          pi_formulas = formula_ispline, mu_formulas = formula_ispline)
# # alpha = 0.1: FDPhat 0.0994, Number of Rej. 332
# # alpha = 0.09: FDPhat 0.0872, Number of Rej. 321
# # alpha = 0.08: FDPhat 0.0784, Number of Rej. 306
#
# # bs
# formula_bs <- c(paste0("bs(log_count, df = ",2:8,")"))
# glm_bs_logcount <- adapt_glm(x, pvals = p_value,  alphas=alphas,
#                     pi_formulas = formula_bs, mu_formulas = formula_bs)
# # alpha = 0.1: FDPhat 0.0997, Number of Rej. 361
# # alpha = 0.09: FDPhat 0.0897, Number of Rej. 301
# # alpha = 0.08: FDPhat 0.0791, Number of Rej. 215

# GMM - p
# ns
gmm_ns_p_logcount <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
                      beta_formulas = formula)

# alpha = 0.1: FDPhat 0.0998, Number of Rej. 411
# alpha = 0.05: FDPhat 0.049, Number of Rej. 286

# # ispline
# gmm_ispline_p_logcount <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
#                            beta_formulas = formula_ispline)
# # alpha = 0.1: FDPhat 0.099, Number of Rej. 414
# # alpha = 0.05: FDPhat 0.0496, Number of Rej. 272
#
# # bs
# gmm_bs_p_logcount <- adapt_gmm(x, pvals = p_value,  alphas=alphas,
#                       beta_formulas = formula_bs)
# # alpha = 0.1: FDPhat 0.099, Number of Rej. 414
# # alpha = 0.05: FDPhat 0.0496, Number of Rej. 282


# GMM - Z
gmm_ns_z_logcount <- adapt_gmm(x, z = z_value, alphas=alphas,
                      beta_formulas = formula, testing = "two_sided")
# alpha = 0.1: FDPhat 0.0998, Number of Rej. 456
# alpha = 0.05: FDPhat 0.0495, Number of Rej. 303

# gmm_bs_z_logcount <- adapt_gmm(x, z = z_value, alphas=alphas,
#                       beta_formulas = formula_bs, testing = "two_sided")
# # alpha = 0.1: FDPhat 0.0989, Number of Rej. 450
# # alpha = 0.05: FDPhat 0.0488, Number of Rej. 297

################## OrderShapeEM #######################

require(OrderShapeEM)
orderfdr_logcount <- OrderShapeEM(pvals = p_value, order.var = log_count,
                         OrderShapeEM.control(trace = TRUE))
sum(orderfdr_logcount$fdr <= 0.05)
# 130

######################## FDRreg ############################
fdr_theo_logcount <- FDRreg(z_value, as.matrix(log_count), nulltype = 'theoretical')
length(which(fdr_theo_logcount$FDR <= 0.05))
# 163

fdr_emp_logcount <- FDRreg(z_value, as.matrix(log_count), nulltype = 'empirical')
length(which(fdr_emp_logcount$FDR <= 0.05))
# 8

######################## IHW ############################
# ihw_01_logcount <- ihw(pvalues = p_value, covariates = log_count, alpha = 0.01)
# ihw_02_logcount <- ihw(pvalues = p_value, covariates = log_count, alpha = 0.02)
# ihw_03_logcount <- ihw(pvalues = p_value, covariates = log_count, alpha = 0.03)
# ihw_04_logcount <- ihw(pvalues = p_value, covariates = log_count, alpha = 0.04)
# ihw_05_logcount <- ihw(pvalues = p_value, covariates = log_count, alpha = 0.05)
# ihw_06_logcount <- ihw(pvalues = p_value, covariates = log_count, alpha = 0.06)
# ihw_07_logcount <- ihw(pvalues = p_value, covariates = log_count, alpha = 0.07)
# ihw_08_logcount <- ihw(pvalues = p_value, covariates = log_count, alpha = 0.08)
# ihw_09_logcount <- ihw(pvalues = p_value, covariates = log_count, alpha = 0.09)
# ihw_1_logcount <- ihw(pvalues = p_value, covariates = log_count, alpha = 0.1)

rejections(ihw_05_logcount)  # 178
rejections(ihw_06_logcount)  # 209
rejections(ihw_07_logcount)  # 247
rejections(ihw_08_logcount)  # 268
rejections(ihw_09_logcount)  # 293
rejections(ihw_1_logcount)   # 316


###########################################################
# Summarize into a different q level table
q_levels <- alphas

# chai
chai_logcount_rejs <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_logcount$lFDR, q, 1))
})

# adapt_glm
adapt_glm_ns_logcount_rej <- glm_ns_logcount$nrejs[1:10]
# adapt_glm_bs_logcount_rej <- glm_bs_logcount$nrejs[1:10]

# adapt_gmm_p
adapt_gmm_p_ns_logcount_rej <- gmm_ns_p_logcount$nrejs[1:10]
# adapt_gmm_p_bs_logcount_rej <- gmm_bs_p_logcount$nrejs[1:10]

# adapt_gmm_z
adapt_gmm_z_ns_logcount_rej <- gmm_ns_z_logcount$nrejs[1:10]
# adapt_gmm_z_bs_logcount_rej <- gmm_bs_z_logcount$nrejs[1:10]

# OrderShapeEM
ordershapeem_logcount_rejs <- sapply(q_levels, function(q) {
  sum(orderfdr_logcount$fdr <= q)
})

# FDRreg
rdrreg_logcount_rejs <- sapply(q_levels, function(q) {
  length(which(fdr_theo_logcount$FDR <= q))
})

# IHW
ihw_logcount_rejs <- c(rejections(ihw_01_logcount), rejections(ihw_02_logcount), rejections(ihw_03_logcount),
              rejections(ihw_04_logcount), rejections(ihw_05_logcount), rejections(ihw_06_logcount),
              rejections(ihw_07_logcount), rejections(ihw_08_logcount), rejections(ihw_09_logcount),
              rejections(ihw_1_logcount))
# [1]  40  69 123 152 178 209 247 268 293 316

# combine into data frame
eraw_sum_logcount <- data.frame(
  q = q_levels,
  chai = chai_logcount_rejs,
  adapt_glm = adapt_glm_ns_logcount_rej,
  # adapt_glm_bs = adapt_glm_bs_logcount_rej,
  adapt_gmm_p = adapt_gmm_p_ns_logcount_rej,
  # adapt_gmm_p_bs = adapt_gmm_p_bs_logcount_rej,
  adapt_gmm_z = adapt_gmm_z_ns_logcount_rej,
  # adapt_gmm_z_bs = adapt_gmm_z_bs_logcount_rej,
  OrderShapeEM = ordershapeem_logcount_rejs,
  FDRreg = rdrreg_logcount_rejs,
  IHW = ihw_logcount_rejs,
  BH = bh_rejs
)




###########################################################
# Plots
library(ggplot2)
library(tidyr)

# into long format
# eraw_long_logcount <- eraw_sum_logcount %>%
#   pivot_longer(
#     cols = -q,
#     names_to = "Method",
#     values_to = "Rejections"
#   )

# into long format
eraw_long_logcount <- eraw_sum_logcount %>%
  pivot_longer(
    cols = -q,
    names_to = "Method",
    values_to = "Rejections"
  ) %>%
  mutate(Method = as.character(Method))
#
# order <- c("chai", "adapt_glm", "adapt_gmm_z", "adapt_gmm_p",
#            "FDRreg", "OrderShapeEM", "IHW", "BH")
#
# eraw_long$Method <- factor(eraw_long$Method, levels = order)

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

p_eraw_logcount <- plot_rejections_vs_q(eraw_long_logcount,
                                   title = "Number of Discoveries vs q, log of normalized counts",
                                   x_breaks = sort(unique(eraw_long_logcount$q)))

pdf("./plots/Erawijantari/Erawijantari_logcount_updated.pdf", width = 8, height = 4)
plot(p_eraw_logcount)
dev.off()

# eraw_chai_logcount <- filter(eraw_long_logcount, Method == "chai")
# eraw_others_logcount <- filter(eraw_long_logcount, Method != "chai")
#
# # plot
# library(ggrepel)
# library(viridisLite)
#
# # color map: viridis for others, red for chai
# methods <- unique(eraw_long_logcount$Method)
# others  <- setdiff(methods, "chai")
# pal     <- viridis(length(others), option = "plasma")
# method_cols <- c(setNames(pal, others), "chai" = "red")
#
# # last point (max q) for each method
# endpoints <- eraw_long_logcount |>
#   dplyr::arrange(Method, q) |>
#   dplyr::group_by(Method) |>
#   dplyr::slice_tail(n = 1) |>
#   dplyr::ungroup()
#
# # plot with end-of-line labels
# p <- ggplot() +
#   geom_line(data = eraw_others_logcount, aes(q, Rejections, color = Method), size = 0.9) +
#   geom_point(data = eraw_others_logcount, aes(q, Rejections, color = Method), size = 1.8) +
#   geom_line(data = eraw_chai_logcount,   aes(q, Rejections, color = Method), size = 1.5) +
#   geom_point(data = eraw_chai_logcount,  aes(q, Rejections, color = Method), size = 2.5) +
#   geom_text_repel(
#     data = endpoints,
#     aes(q, Rejections, label = Method, color = Method),
#     nudge_x = diff(range(eraw_long_logcount$q)) * 0.03,  # push labels a bit to the right
#     hjust = 0, direction = "y",
#     segment.size = 0.2, box.padding = 0.1, point.padding = 0.1,
#     show.legend = FALSE
#   ) +
#   scale_color_manual(values = method_cols) +
#   labs(title = "Number of Discoveries vs q", x = "q level", y = "Number of Rejections") +
#   coord_cartesian(clip = "off") +
#   expand_limits(x = max(eraw_long_logcount$q) + 0.01) +         # room for labels
#   theme_minimal() +
#   theme(
#     plot.margin = margin(5.5, 60, 5.5, 5.5),           # extra right margin
#     legend.position = "none"                            # labels replace legend; set to "bottom" if you want both
#   )

# save(eraw_sum, eraw_sum_logcount, file = "./chai_env/Erawijantari_tables.RData", compress = "xz")




