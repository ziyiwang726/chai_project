# AE and PCA are from python ONLY

# ZhuF_2020
suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(GGally)
})

# ------------------------------
# Read autoencoder and PCA5 file
# ------------------------------
ae5  <- read.csv("./data/Zhu_autoencoder/outputs/ZhuF_2020_species_AE5.csv", row.names = 1)
pca5 <- read.csv("./data/Zhu_autoencoder/outputs/ZhuF_2020_species_PCA5.csv", row.names = 1)
# spec_path <- read.csv("./data/Zhu_autoencoder/ZhuF_2020_species_pathway.csv", row.names = 1)


# -----------------------------
# Helper: extract species label
# -----------------------------
.extract_species <- function(taxon) {
  sp <- ifelse(grepl("s__", taxon, fixed = TRUE),
               sub(".*s__([^|;]+).*", "\\1", taxon),
               taxon)
  sp <- gsub("^s__", "", sp)
  sp <- gsub("_", " ", sp)
  trimws(sp)
}


# ------------------------------------------------
# Extract relative abundance table
# Apply the same rule of relative abundance table
# ------------------------------------------------
# Extract the relative abundance
abun <- curatedMetagenomicData("ZhuF_2020.relative_abundance", dryrun = FALSE)

outcome <- as.factor(abun$`2021-03-31.ZhuF_2020.relative_abundance`[[5]])
# schizophrenia 90 vs. control 81

# Extract the assay matrix
abun_mat <- assay(abun[[1]], "relative_abundance")
# 480 x 171

r_names <- rownames(abun_mat)
taxa_names <- sub("^.*\\|", "", r_names)
species_clean <- .extract_species(taxa_names)
rownames(abun_mat) <- species_clean
# 480 x 171

# Check whether subject order is the same
table(abun$`2021-03-31.ZhuF_2020.relative_abundance`[[2]] == colnames(abun_mat))


# -------------------------
# Keep sharing species only
# -------------------------

length(intersect(rownames(ae5), rownames(abun_mat)))
common_species <- intersect(rownames(ae5), rownames(abun_mat))
# 323 sharing

abun_mat_common <- abun_mat[common_species, , drop = FALSE]

# Column-Normalize (column is sample)
cs <- colSums(abun_mat_common, na.rm = TRUE)
cs[cs == 0] <- 1
abun_mat_norm <- sweep(abun_mat_common, 2, cs, FUN = "/")

t_abun <- t(abun_mat_norm)
# 171 subjects x 323 species


# Checking the order of species
table(colnames(t_abun) == rownames(ae5))
# 323 species x 5 AEs

table(colnames(t_abun) == rownames(pca5))
# 323 species x 5 PCs

# ----------------
# p-value/z-value
# ----------------
library(coin)
p_z_values <- sapply(colnames(t_abun), function(feature) {
  features <- t_abun[, feature, drop = FALSE]
  wilcox <- wilcox_test(features[,1] ~ outcome)

  p_value <- coin::pvalue(wilcox)
  z_value <- statistic(wilcox)

  return(c(p_value = p_value, z_value = z_value))
}, simplify = "matrix")
p_z_values <- t(p_z_values)

p_value <- p_z_values[,1]
z_value <- p_z_values[,2]

##################### BH ###########################
p_adjusted <- p.adjust(p_value, method = "BH")
bh <- which(p_adjusted <= 0.05)
length(bh) # 4


##################### Visualizaiton ######################
library(ggplot2)
library(viridis)

df <- data.frame(z = z_value, ae5)
# Latent space plot
ggplot(df, aes(x = AE1, y = AE2, color = z)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_viridis(option = "D") +
  theme_minimal() +
  labs(
    title = "Autoencoder Latent Space Map",
    subtitle = "Position defined by AE dimensions AE1 and AE2; color defined by z-statistics",
    x = "Autoencoder Dimension 1",
    y = "Autoencoder Dimension 2",
    color = "z_value"
  ) +
  theme(legend.position = "right")



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

source("chai.R")
source("color_helper.R")

library(IHW)
library(FDRreg)


################## random X #######################
set.seed(123)
# Standard normal noise: N(0, 1)
noise <- rnorm(length(z_value))

# # Uniform noise on [0, 1]
noise2 <- runif(length(z_value))

################## chai #######################
set.seed(123)
chai_zhu_ae5 <- chai(z_value, ae5, M = 100)

length(lFDRselect(chai_zhu_ae5$lFDR, 0.01, 1)) # 6
length(lFDRselect(chai_zhu_ae5$lFDR, 0.05, 1)) # 38
length(lFDRselect(chai_zhu_ae5$lFDR, 0.1, 1))  # 65


chai_zhu_pca5 <- chai(z_value, pca5, M = 100)

length(lFDRselect(chai_zhu_pca5$lFDR, 0.01, 1)) # 6
length(lFDRselect(chai_zhu_pca5$lFDR, 0.05, 1)) # 29
length(lFDRselect(chai_zhu_pca5$lFDR, 0.1, 1))  # 55

chai_zhu_noise <- chai(z_value, noise, M = 100)
length(lFDRselect(chai_zhu_noise$lFDR, 0.01, 1)) # 3
length(lFDRselect(chai_zhu_noise$lFDR, 0.05, 1)) # 20
length(lFDRselect(chai_zhu_noise$lFDR, 0.1, 1))  # 53

chai_zhu_noise2 <- chai(z_value, noise2, M = 100)


length(intersect(which(p_adjusted < 0.05), lFDRselect(chai_zhu_ae5$lFDR, 0.05, 1)))
length(intersect(which(p_adjusted < 0.1), lFDRselect(chai_zhu_ae5$lFDR, 0.1, 1)))


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

chai_sel_ae5 <- label_direction(outcome, z_value, lFDRselect(chai_zhu_ae5$lFDR, 0.05, 1))
# chai_sel$feature[chai_sel$direction == "higher_in_cirrhosis"]


# Visualization
chai_sel <- lFDRselect(chai_zhu_ae5$lFDR, 0.05, 1)
df2 <- data.frame(z_value, ae5) %>%
  dplyr::mutate(
    grp = dplyr::case_when(
      dplyr::row_number() %in% bh & dplyr::row_number() %in% chai_sel ~ "both",
      dplyr::row_number() %in% bh ~ "BH_only",
      dplyr::row_number() %in% chai_sel ~ "chai_only",
      TRUE ~ "neither"
    )
  )

ggplot(df2, aes(x = AE1, y = AE2, color = grp)) +
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
    title = "Autoencoder Latent Space Map",
    subtitle = "Points colored by findings",
    x = "Autoencoder Dimension 1",
    y = "Autoencoder Dimension 2",
    color = "Group"
  ) +
  theme(legend.position = "right")

library(VennDiagram)
ae5_sel <- lFDRselect(chai_zhu_ae5$lFDR, 0.05, 1)
pca5_sel <- lFDRselect(chai_zhu_pca5$lFDR, 0.05, 1)
noise1_sel <- lFDRselect(chai_zhu_noise$lFDR, 0.05, 1)
noise2_sel <- lFDRselect(chai_zhu_noise2$lFDR, 0.05, 1)
# bh

# Chart
venn.diagram(
  x = list(ae5_sel, pca5_sel, noise1_sel, noise2_sel, bh),
  category.names = c("Autoencoder" , "PCA" , "noise 1", "noise 2", "BH"),
  filename = 'ae_pca_venn_diagramm.png',
  output=TRUE
)

# UpSet
ae5_sel_tax <- names(z_value)[ae5_sel]
pca5_sel_tax <- names(z_value)[pca5_sel]
noise1_sel_tax <- names(z_value)[noise1_sel]
noise2_sel_tax <- names(z_value)[noise2_sel]
bh_sel_tax <- names(z_value)[bh]

all_species <- sort(unique(c(ae5_sel_tax, pca5_sel_tax, noise1_sel_tax, noise2_sel_tax, bh_sel_tax)))

upset_df <- tibble(species = all_species) %>%
  mutate(
    AE     = species %in% ae5_sel_tax,
    PCA    = species %in% pca5_sel_tax,
    noise1 = species %in% noise1_sel_tax,
    noise2 = species %in% noise2_sel_tax,
    BH     = species %in% bh_sel_tax
  )

upset(upset_df,
      intersect = c("AE","PCA","noise1","noise2","BH"),
      base_annotations = list(
        "Intersection size" = intersection_size())
      )

################## Adapt #######################
# alphas <- c(0.001, 0.01, 0.05)
alphas <- seq(0.01, 0.1, by = 0.01)


# GLM
# ns - AE5
x_ae5 <- data.frame(ae5)
make_ns_formulas <- function(dfs) {
  lapply(dfs, function(d)
    as.formula(paste0("~ ns(AE1, df=", d, ") + ns(AE2, df=", d, ") + ns(AE3, df=", d, ") +
                      ns(AE4, df=", d, ") + ns(AE5, df=", d, ")"))
  )
}
dfs <- c(2, 4, 6, 8)
formula_ns_ae5 <- make_ns_formulas(dfs)

glm_ns_ae5 <- adapt_glm(x_ae5, pvals = p_value, alphas=alphas,
                        pi_formulas = formula_ns_ae5, mu_formulas = formula_ns_ae5)
# alpha = 0.1: FDPhat 0.0882, Number of Rej. 34
# alpha = 0.09: FDPhat 0.0882, Number of Rej. 34
# alpha = 0.08: FDPhat 0.069, Number of Rej. 29
# alpha = 0.07: FDPhat 0.069, Number of Rej. 29

length(intersect(which(p_adjusted <= 0.1), glm_ns_ae5$rejs[10][[1]]))

# ns - PCA5
x_pca5 <- data.frame(pca5)
make_ns_formulas <- function(dfs) {
  lapply(dfs, function(d)
    as.formula(paste0("~ ns(PC1, df=", d, ") + ns(PC2, df=", d, ") + ns(PC3, df=", d, ") +
                      ns(PC4, df=", d, ") + ns(PC5, df=", d, ")"))
  )
}
dfs <- c(2, 4, 6, 8)
formula_ns_pca5 <- make_ns_formulas(dfs)

glm_ns_pca5 <- adapt_glm(x_pca5, pvals = p_value, alphas=alphas,
                         pi_formulas = formula_ns_pca5, mu_formulas = formula_ns_pca5)
# alpha = 0.1: FDPhat 0.087, Number of Rej. 46
# alpha = 0.09: FDPhat 0.087, Number of Rej. 46
# alpha = 0.08: FDPhat 0.0698, Number of Rej. 43
# alpha = 0.07: FDPhat 0.0698, Number of Rej. 43
# alpha = 0.06: FDPhat 0.0541, Number of Rej. 37

length(intersect(which(p_adjusted <= 0.1), glm_ns_pca5$rejs[10][[1]]))

# ns - noise
x_noise <- data.frame(noise)
formula_ns_noise <- paste0("ns(noise, df = ", 2:8, ")")

glm_ns_noise <- adapt_glm(x_noise, pvals = p_value, alphas=alphas,
                          pi_formulas = formula_ns_noise, mu_formulas = formula_ns_noise)
# alpha = 0.1: FDPhat 0.08, Number of Rej. 50
# alpha = 0.09: FDPhat 0.08, Number of Rej. 50
# alpha = 0.08: FDPhat 0.08, Number of Rej. 50
# alpha = 0.07: FDPhat 0.06, Number of Rej. 50
# alpha = 0.06: FDPhat 0.06, Number of Rej. 50

# ns - noise2
x_noise2 <- data.frame(noise2)
formula_ns_noise2 <- paste0("ns(noise2, df = ", 2:8, ")")

glm_ns_noise2 <- adapt_glm(x_noise2, pvals = p_value, alphas=alphas,
                           pi_formulas = formula_ns_noise2, mu_formulas = formula_ns_noise2)
# alpha = 0.1: FDPhat 0.0943, Number of Rej. 53
# alpha = 0.09: FDPhat 0.0769, Number of Rej. 52
# alpha = 0.08: FDPhat 0.0769, Number of Rej. 52
# alpha = 0.07: FDPhat 0.0588, Number of Rej. 51
# alpha = 0.06: FDPhat 0.0588, Number of Rej. 51


# GMM - p
# ns - AE5
gmm_ns_p_ae5 <- adapt_gmm(x_ae5, pvals = p_value,  alphas=alphas,
                          beta_formulas = formula_ns_ae5)
# alpha = 0.1: FDPhat 0.0997, Number of Rej. 54
# alpha = 0.05: FDPhat 0.0449, Number of Rej. 6
# alpha = 0.01: FDPhat 0.0179, Number of Rej. 0

length(intersect(which(p_adjusted <= 0.05), gmm_ns_p_ae5$rejs[5][[1]]))

# ns - PCA5
gmm_ns_p_pca5 <- adapt_gmm(x_pca5, pvals = p_value,  alphas=alphas,
                           beta_formulas = formula_ns_pca5)
# alpha = 0.1: FDPhat 0.0997, Number of Rej. 54
# alpha = 0.05: FDPhat 0.0449, Number of Rej. 6
# alpha = 0.01: FDPhat 0.0179, Number of Rej. 0

length(intersect(which(p_adjusted <= 0.05), gmm_ns_p_pca5$rejs[5][[1]]))

# ns - noise
gmm_ns_p_noise <- adapt_gmm(x_noise, pvals = p_value,  alphas=alphas,
                            beta_formulas = formula_ns_noise)
# alpha = 0.1: FDPhat 0.0995, Number of Rej. 53
# alpha = 0.05: FDPhat 0.0449, Number of Rej. 6
# alpha = 0.01: FDPhat 0.0179, Number of Rej. 0

# ns - noise2
gmm_ns_p_noise2 <- adapt_gmm(x_noise2, pvals = p_value,  alphas=alphas,
                             beta_formulas = formula_ns_noise2)
# alpha = 0.1: FDPhat 0.0995, Number of Rej. 53
# alpha = 0.05: FDPhat 0.0449, Number of Rej. 6
# alpha = 0.01: FDPhat 0.0179, Number of Rej. 0

# GMM - Z - AE5
gmm_ns_z_ae5 <- adapt_gmm(x_ae5, z = z_value, alphas=alphas,
                          beta_formulas = formula_ns_ae5, testing = "two_sided")
# alpha = 0.1: FDPhat 0.1, Number of Rej. 7
# alpha = 0.05: FDPhat 0.0269, Number of Rej. 2
# alpha = 0.01: FDPhat 0.0269, Number of Rej. 0

length(intersect(which(p_adjusted <= 0.05), gmm_ns_z_ae5$rejs[5][[1]]))


# GMM - Z - PCA5
gmm_ns_z_pca5 <- adapt_gmm(x_pca5, z = z_value, alphas=alphas,
                           beta_formulas = formula_ns_pca5, testing = "two_sided")
# alpha = 0.1: FDPhat 0.1, Number of Rej. 7
# alpha = 0.05: FDPhat 0.0269, Number of Rej. 2
# alpha = 0.01: FDPhat 0.0269, Number of Rej. 0

length(intersect(which(p_adjusted <= 0.05), gmm_ns_z_pca5$rejs[5][[1]]))

# ns - noise
gmm_ns_z_noise <- adapt_gmm(x_noise, z = z_value,  alphas=alphas,
                            beta_formulas = formula_ns_noise, testing = "two_sided")
# alpha = 0.1: FDPhat 0.0969, Number of Rej. 10
# alpha = 0.05: FDPhat 0.0359, Number of Rej. 3
# alpha = 0.01: FDPhat 0.0179, Number of Rej. 0

# ns - noise
gmm_ns_z_noise2 <- adapt_gmm(x_noise2, z = z_value,  alphas=alphas,
                             beta_formulas = formula_ns_noise2, testing = "two_sided")
# alpha = 0.1: FDPhat 0.1, Number of Rej. 7
# alpha = 0.05: FDPhat 0.0269, Number of Rej. 2
# alpha = 0.01: FDPhat 0.0179, Number of Rej. 0

################## OrderShapeEM #######################
require(OrderShapeEM)
ordershape_x1_ae5 <- OrderShapeEM(pvals = p_value, order.var = ae5[,1],
                                  OrderShapeEM.control(trace = TRUE))
# sum(ordershape_x1$fdr <= 0.001) # 2
sum(ordershape_x1_ae5$fdr <= 0.01)  # 4
sum(ordershape_x1_ae5$fdr <= 0.05)  # 30
sum(ordershape_x1_ae5$fdr <= 0.1)   # 64

length(intersect(which(p_adjusted <= alphas[5]),  which(ordershape_x1_ae5$fdr <= alphas[5])))


ordershape_noise <- OrderShapeEM(pvals = p_value, order.var = noise,
                                 OrderShapeEM.control(trace = TRUE))
sum(ordershape_noise$fdr <= 0.01)  # 2
sum(ordershape_noise$fdr <= 0.05)  # 14
sum(ordershape_noise$fdr <= 0.1)   # 53

ordershape_noise2 <- OrderShapeEM(pvals = p_value, order.var = noise2,
                                  OrderShapeEM.control(trace = TRUE))
sum(ordershape_noise2$fdr <= 0.01)  # 2
sum(ordershape_noise2$fdr <= 0.05)  # 14
sum(ordershape_noise2$fdr <= 0.1)   # 53


# mat: ae5_common or pca5_common
run_ordershape_by_col <- function(mat, dataset = c("ae5", "pca5")) {
  dataset <- match.arg(dataset)

  res <- vector("list", ncol(mat))
  names(res) <- paste0("ordershape_x", seq_len(ncol(mat)), "_", dataset)

  for (j in seq_len(ncol(mat))) {
    message("Fitting OrderShapeEM for column ", j)
    res[[j]] <- OrderShapeEM(
      pvals     = p_value,
      order.var = mat[, j],
      OrderShapeEM.control(trace = TRUE)
    )
  }
  res
}

## For ae5_common (5 columns)
ordershape_ae5 <- run_ordershape_by_col(ae5, dataset = "ae5")

## For ae5_common (5 columns)
ordershape_pca5 <- run_ordershape_by_col(pca5, dataset = "pca5")



######################## FDRreg ############################
set.seed(123)
fdr_theo_ae5 <- FDRreg(z_value, as.matrix(ae5), nulltype = 'theoretical')
length(which(fdr_theo_ae5$FDR <= 0.01))   # 2
length(which(fdr_theo_ae5$FDR <= 0.05))   # 25
length(which(fdr_theo_ae5$FDR <= 0.1))    # 58

length(intersect(which(p_adjusted <= alphas[5]),  which(fdr_theo_ae5$FDR <= alphas[5])))

fdr_theo_pca5 <- FDRreg(z_value, as.matrix(pca5), nulltype = 'theoretical')
length(which(fdr_theo_pca5$FDR <= 0.01))   # 3
length(which(fdr_theo_pca5$FDR <= 0.05))   # 23
length(which(fdr_theo_pca5$FDR <= 0.1))    # 58


fdr_theo_noise <- FDRreg(z_value, as.matrix(noise), nulltype = 'theoretical')
length(which(fdr_theo_noise$FDR <= 0.01))   # 2
length(which(fdr_theo_noise$FDR <= 0.05))   # 19
length(which(fdr_theo_noise$FDR <= 0.1))    # 55

fdr_theo_noise2 <- FDRreg(z_value, as.matrix(noise2), nulltype = 'theoretical')
length(which(fdr_theo_noise2$FDR <= 0.01))   # 2
length(which(fdr_theo_noise2$FDR <= 0.05))   # 18
length(which(fdr_theo_noise2$FDR <= 0.1))    # 55

###########################################################
# Summarize into a different q level table
q_levels <- alphas

# chai
chai_rejs_ae5 <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_zhu_ae5$lFDR, q, 1))
})

chai_rejs_pca5 <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_zhu_pca5$lFDR, q, 1))
})

chai_rejs_noise <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_zhu_noise$lFDR, q, 1))
})

chai_rejs_noise2 <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_zhu_noise2$lFDR, q, 1))
})

# adapt_glm
adapt_glm_ns_ae5 <- glm_ns_ae5$nrejs
adapt_glm_ns_pca5 <- glm_ns_pca5$nrejs
adapt_glm_ns_noise <- glm_ns_noise$nrejs
adapt_glm_ns_noise2 <- glm_ns_noise2$nrejs


# adapt_gmm_p
adapt_gmm_p_ns_ae5 <- gmm_ns_p_ae5$nrejs
adapt_gmm_p_ns_pca5 <- gmm_ns_p_pca5$nrejs
adapt_gmm_p_ns_noise <- gmm_ns_p_noise$nrejs
adapt_gmm_p_ns_noise2 <- gmm_ns_p_noise2$nrejs


# adapt_gmm_z
adapt_gmm_z_ns_ae5 <- gmm_ns_z_ae5$nrejs
adapt_gmm_z_ns_pca5 <- gmm_ns_z_pca5$nrejs
adapt_gmm_z_ns_noise <- gmm_ns_z_noise$nrejs
adapt_gmm_z_ns_noise2 <- gmm_ns_z_noise2$nrejs


# OrderShapeEM
count_ordershape_selections <- function(os_list, q_levels) {
  # os_list: list like ordershape_ae5 or ordershape_pca5
  # returns a matrix: rows = q_levels, cols = each x1/x2/...
  res <- sapply(os_list, function(os) {
    # for one ordershape object
    sapply(q_levels, function(q) sum(os$fdr <= q))
  })
  # sapply gives q_levels in rows; make that explicit:
  rownames(res) <- paste0("q=", q_levels)
  res
}

ordershape_ae5  <- count_ordershape_selections(ordershape_ae5,  q_levels)
ordershape_pca5 <- count_ordershape_selections(ordershape_pca5, q_levels)

ordershapeem_noise <- sapply(q_levels, function(q) {
  sum(ordershape_noise$fdr <= q)
})

ordershapeem_noise2 <- sapply(q_levels, function(q) {
  sum(ordershape_noise2$fdr <= q)
})

# FDRreg
rdrreg_ae5 <- sapply(q_levels, function(q) {
  length(which(fdr_theo_ae5$FDR <= q))
})

rdrreg_pca5 <- sapply(q_levels, function(q) {
  length(which(fdr_theo_pca5$FDR <= q))
})

rdrreg_noise <- sapply(q_levels, function(q) {
  length(which(fdr_theo_noise$FDR <= q))
})

rdrreg_noise2 <- sapply(q_levels, function(q) {
  length(which(fdr_theo_noise2$FDR <= q))
})

# IHW
# run_ihw_by_col <- function(mat, q_levels, dataset = c("ae5", "pca5"), nbins = 5) {
#   dataset <- match.arg(dataset)
#
#   # result: rows = q_levels, cols = columns of mat
#   res_mat <- sapply(seq_len(ncol(mat)), function(j) {
#     sapply(q_levels, function(q) {
#       fit <- ihw(
#         pvalues   = p_value,
#         covariates = mat[, j],
#         alpha     = q,
#         nbins     = nbins
#       )
#       rejections(fit)
#     })
#   })
#
#   # sapply puts q_levels in rows; make that explicit
#   rownames(res_mat) <- paste0("q=", q_levels)
#   colnames(res_mat) <- paste0("ihw_x", seq_len(ncol(mat)), "_", dataset)
#
#   res_mat
# }
#
# ihw_ae5  <- run_ihw_by_col(ae5_common,  q_levels, dataset = "ae5")
# ihw_pca5 <- run_ihw_by_col(pca5_common, q_levels, dataset = "pca5")

# Careful on IHW!!!!! May cause R terminated!
ihw_x1_ae5 <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = ae5[,1], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_x2_ae5 <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = ae5[,2], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_x3_ae5 <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = ae5[,3], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_x4_ae5 <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = ae5[,4], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_x5_ae5 <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = ae5[,5], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_x1_pca5 <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = pca5[,1], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_x2_pca5 <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = pca5[,2], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_x3_pca5 <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = pca5[,3], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_x4_pca5 <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = pca5[,4], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_x5_pca5 <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = pca5[,5], alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_noise <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = noise, alpha = q, nbins = 5)
  rejections(ihw)
})

ihw_noise2 <- sapply(q_levels, function(q) {
  ihw <- ihw(pvalues = p_value, covariates = noise2, alpha = q, nbins = 5)
  rejections(ihw)
})

# BH
bh_rejs <- sapply(q_levels, function(q) {
  length(which(p_adjusted <= q))
})


# combine into data frame
zhu_sum_ae5 <- data.frame(
  q = q_levels,
  chai = chai_rejs_ae5,
  adapt_glm = adapt_glm_ns_ae5,
  adapt_gmm_p = adapt_gmm_p_ns_ae5,
  adapt_gmm_z = adapt_gmm_z_ns_ae5,
  OrderShapeEM_x1 = ordershape_ae5[, 1],
  OrderShapeEM_x2 = ordershape_ae5[, 2],
  OrderShapeEM_x3 = ordershape_ae5[, 3],
  OrderShapeEM_x4 = ordershape_ae5[, 4],
  OrderShapeEM_x5 = ordershape_ae5[, 5],
  FDRreg = rdrreg_ae5,
  IHW_x1 = ihw_x1_ae5,
  IHW_x2 = ihw_x2_ae5,
  IHW_x3 = ihw_x3_ae5,
  IHW_x4 = ihw_x4_ae5,
  IHW_x5 = ihw_x5_ae5,
  BH = bh_rejs
)

zhu_sum_pca5 <- data.frame(
  q = q_levels,
  chai = chai_rejs_pca5,
  adapt_glm = adapt_glm_ns_pca5,
  adapt_gmm_p = adapt_gmm_p_ns_pca5,
  adapt_gmm_z = adapt_gmm_z_ns_pca5,
  OrderShapeEM_x1 = ordershape_pca5[, 1],
  OrderShapeEM_x2 = ordershape_pca5[, 2],
  OrderShapeEM_x3 = ordershape_pca5[, 3],
  OrderShapeEM_x4 = ordershape_pca5[, 4],
  OrderShapeEM_x5 = ordershape_pca5[, 5],
  FDRreg = rdrreg_pca5,
  IHW_x1 = ihw_x1_pca5,
  IHW_x2 = ihw_x2_pca5,
  IHW_x3 = ihw_x3_pca5,
  IHW_x4 = ihw_x4_pca5,
  IHW_x5 = ihw_x5_pca5,
  BH = bh_rejs
)


###########################################################
# Plots
library(ggplot2)
library(tidyr)
library(ggrepel)
library(viridisLite)

# AE
zhu_ae <- zhu_sum_ae5 %>%
  pivot_longer(
    cols = -q,
    names_to = "Method",
    values_to = "Rejections"
  ) %>%
  mutate(Method = as.character(Method))

p_zhu_ae <- plot_rejections_vs_q(zhu_ae,
                              title = "Number of Discoveries vs q, with Autoencoder",
                              x_breaks = sort(unique(zhu_ae$q)))

# pdf("./plots/Zhu/Zhu_AE_updated.pdf", width = 10, height = 5)
# plot(p_zhu_ae)
# dev.off()


# PCA
zhu_pca <- zhu_sum_pca5 %>%
  pivot_longer(
    cols = -q,
    names_to = "Method",
    values_to = "Rejections"
  ) %>%
  mutate(Method = as.character(Method))

p_zhu_pca <- plot_rejections_vs_q(zhu_pca,
                                 title = "Number of Discoveries vs q, with PCA",
                                 x_breaks = sort(unique(zhu_pca$q)))

pdf("./plots/Zhu/Zhu_PCA_updated.pdf", width = 10, height = 5)
plot(p_zhu_pca)
dev.off()



######################################### AE vs. PCA #################################################
zhu_sum_ae_pca <- data.frame(
  q = q_levels,
  chai_ae5 = chai_rejs_ae5,
  chai_pca5 = chai_rejs_pca5,
  chai_noise1 = chai_rejs_noise,
  chai_noise2 = chai_rejs_noise2,
  adapt_glm_ae5 = adapt_glm_ns_ae5,
  adapt_glm_pca5 = adapt_glm_ns_pca5,
  adapt_glm_noise1 = adapt_glm_ns_noise,
  adapt_glm_noise2 = adapt_glm_ns_noise2,
  adapt_gmm_p_ae5 = adapt_gmm_p_ns_ae5,
  adapt_gmm_p_pca5 = adapt_gmm_p_ns_pca5,
  adapt_gmm_p_noise1 = adapt_gmm_p_ns_noise,
  adapt_gmm_p_noise2 = adapt_gmm_p_ns_noise2,
  adapt_gmm_z_ae5 = adapt_gmm_z_ns_ae5,
  adapt_gmm_z_pca5 = adapt_gmm_z_ns_pca5,
  adapt_gmm_z_noise1 = adapt_gmm_z_ns_noise,
  adapt_gmm_z_noise2 = adapt_gmm_z_ns_noise2,
  fdrreg_ae5 = rdrreg_ae5,
  fdrreg_pca5 = rdrreg_pca5,
  fdrreg_noise1 = rdrreg_noise,
  fdrreg_noise2 = rdrreg_noise2,
  bh = bh_rejs
)

###########################################################
# Plots
zhu_long <- zhu_sum_ae_pca %>%
  pivot_longer(
    cols = -q,
    names_to  = "name",
    values_to = "Rejections"
  ) %>%
  mutate(
    side = case_when(
      grepl("_ae5$",   name) ~ "ae5",
      grepl("_pca5$",  name) ~ "pca5",
      grepl("_noise1$", name) ~ "noise1",
      grepl("_noise2$", name) ~ "noise2",
      TRUE                    ~ "bh"        # BH has no suffix
    ),
    Method = case_when(
      side == "bh" ~ "bh",
      TRUE         ~ sub("_(ae5|pca5|noise1|noise2)$", "", name)
    )
  )

zhu_long$side <- factor(zhu_long$side, levels = c("ae5", "pca5", "noise1", "noise2"))

# Only methods with the three side-info types (drop BH panel)
zhu_long_methods <- zhu_long %>%
  filter(side != "bh", Method != "bh")

# --- Put 'chai' as the first facet ----------------------------------------
method_levels <- c("chai", setdiff(unique(as.character(zhu_long_methods$Method)), "chai"))

zhu_long_methods <- zhu_long_methods %>%
  mutate(Method = factor(Method, levels = method_levels))

# --- End-of-line label positions (one per Method × side) -------------------
endpoints <- zhu_long_methods %>%
  arrange(Method, side, q) %>%
  group_by(Method, side) %>%
  slice_tail(n = 1) %>%
  ungroup()

dx <- diff(range(zhu_long_methods$q))

# BH baseline (as before)
bh_line <- zhu_long %>%
  dplyr::filter(Method == "bh") %>%
  dplyr::select(q, Rejections) %>%
  tidyr::crossing(Method = method_levels) %>%
  # IMPORTANT: make Method a factor with same levels
  mutate(Method = factor(Method, levels = method_levels))

# Endpoints for BH (one per Method)
bh_endpoints <- bh_line %>%
  group_by(Method) %>%
  slice_max(q, with_ties = FALSE) %>%
  ungroup()

p_ae_pca <- ggplot(zhu_long_methods,
                   aes(x = q, y = Rejections, color = side)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_text_repel(
    data = endpoints,
    aes(label = side),
    nudge_x      = dx * 0.03,
    hjust        = 0,
    direction    = "y",
    segment.size = 0.2,
    box.padding  = 0.1,
    point.padding= 0.1,
    show.legend  = FALSE
  ) +
  # BH dashed line
  geom_line(
    data        = bh_line,
    aes(x = q, y = Rejections, group = 1),
    inherit.aes = FALSE,
    color       = "black",
    linewidth   = 0.8,
    linetype    = "dashed"
  ) +
  # "bh" label
  geom_text(
    data        = bh_endpoints,
    aes(x = q, y = Rejections, label = "bh"),
    inherit.aes = FALSE,
    nudge_x     = dx * 0.03,
    hjust       = 0,
    vjust       = 0.5,
    size        = 3,
    color       = "black"
  ) +
  facet_wrap(~ Method) +  # <-- remove scales = "free_y"
  scale_x_continuous(breaks = q_levels) +
  scale_y_continuous(limits = c(0, max(zhu_long_methods$Rejections))) +
  scale_color_manual(
    values = c(
      ae5    = "red",
      pca5   = "blue",
      noise1 = "grey50",
      noise2 = "grey70"
    )
  ) +
  labs(
    title  = "Number of Discoveries at different q level",
    x      = "q level",
    y      = "Number of Rejections",
    color  = "Side information: "
  ) +
  coord_cartesian(clip = "off") +
  expand_limits(x = max(zhu_long_methods$q) + dx * 0.05) +
  theme_minimal() +
  theme(
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin     = margin(5.5, 60, 5.5, 5.5)
  )


# pdf("./plots/Zhu/Zhu_AE5_vs_PCA5_vs_noise.pdf",
#     width = 11, height = 7)
# plot(p_ae_pca)
# dev.off()
