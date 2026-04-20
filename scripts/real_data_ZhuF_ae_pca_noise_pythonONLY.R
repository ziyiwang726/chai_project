# AE and PCA are from python ONLY

# ZhuF_2020
project_lib <- "./r_libs"
if (dir.exists(project_lib)) {
  .libPaths(c(normalizePath(project_lib), .libPaths()))
}

experimenthub_cache <- "./data/ExperimentHub"
dir.create(experimenthub_cache, recursive = TRUE, showWarnings = FALSE)
if (!nzchar(Sys.getenv("EXPERIMENT_HUB_CACHE"))) {
  Sys.setenv(EXPERIMENT_HUB_CACHE = normalizePath(experimenthub_cache, mustWork = FALSE))
}
if (requireNamespace("ExperimentHub", quietly = TRUE)) {
  ExperimentHub::setExperimentHubOption("LOCAL", TRUE)
}

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

dir.create("./plots/Zhu", recursive = TRUE, showWarnings = FALSE)

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
abun <- curatedMetagenomicData(
  "ZhuF_2020.relative_abundance",
  dryrun = FALSE
)

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
library(mclust)
library(adaptMT)
library(AdaPTGMM)
library(splines)

source_first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (!length(existing)) {
    stop("Could not locate required helper file. Tried: ", paste(paths, collapse = ", "))
  }
  source(existing[[1]])
}

source_first_existing(c(
  "./code/chai.R",
  "./Gopalakrishnan/code/chai.R",
  "./Gopalakrishnan03252026/codeZW/chai.R"
))
source_first_existing(c(
  "./code/conditionalParam.R",
  "./Gopalakrishnan/code/conditionalParam.R",
  "./Gopalakrishnan03252026/codeZW/conditionalParam.R"
))
source_first_existing(c(
  "./code/naiveRemoveOneObs.R",
  "./Gopalakrishnan/code/naiveRemoveOneObs.R",
  "./Gopalakrishnan03252026/codeZW/naiveRemoveOneObs.R"
))
source_first_existing(c(
  "./code/rGaussianMix.R",
  "./Gopalakrishnan/code/rGaussianMix.R",
  "./Gopalakrishnan03252026/codeZW/rGaussianMix.R"
))
source_first_existing(c(
  "./code/utils.R",
  "./Gopalakrishnan/code/utils.R",
  "./Gopalakrishnan03252026/codeZW/utils.R"
))
source_first_existing(c(
  "./code/color_helper.R",
  "./Gopalakrishnan/code/color_helper.R",
  "./Gopalakrishnan03252026/codeZW/color_helper.R"
))

# Extend the shared plotting helpers with Efron's theoretical-null locfdr baseline.
.base_make_method_palette <- make_method_palette
.base_method_order <- method_order
efron_method_labels <- c("Efron_lFDR", "Efron lFDR N(0,1)")

make_method_palette <- function(methods, K = 12) {
  methods <- unique(as.character(methods))
  pal <- .base_make_method_palette(methods, K = K)

  present_efron <- intersect(efron_method_labels, methods)
  if (length(present_efron)) {
    pal[present_efron] <- "#4D4D4D"
  }

  pal[methods]
}

method_order <- function(methods) {
  methods <- unique(as.character(methods))
  ord <- .base_method_order(setdiff(methods, efron_method_labels))

  present_efron <- intersect(efron_method_labels, methods)
  if (length(present_efron)) {
    efron_method <- present_efron[[1]]
    bh_idx <- match("BH", ord, nomatch = 0L)
    if (bh_idx > 0L) {
      ord <- append(ord, efron_method, after = bh_idx)
    } else {
      ord <- c(ord, efron_method)
    }
  }

  unique(c(ord, setdiff(methods, ord)))
}

library(IHW)
fdrreg_available <- requireNamespace("FDRreg", quietly = TRUE)
if (!fdrreg_available) {
  message("FDRreg not available; skipping FDRreg fits and summaries.")
}

################## side-info-free baselines #######################
# Use Efron's locfdr with the theoretical N(0,1) null.
efron_fit <- locfdr::locfdr(z_value, plot = 0, nulltype = 0)
efron_lfdr <- as.numeric(efron_fit$fdr)
efron_lfdr[!is.finite(efron_lfdr)] <- 1
efron_lfdr <- pmin(pmax(efron_lfdr, 0), 1)
names(efron_lfdr) <- names(z_value)

efron <- which(efron_lfdr <= 0.05)


################## random X #######################
set.seed(123)
# Standard normal noise: N(0, 1)
noise <- rnorm(length(z_value))

# # Uniform noise on [0, 1]
noise2 <- runif(length(z_value))

# Five-dimensional Gaussian noise side information
noise5 <- matrix(rnorm(length(z_value) * 5), nrow = length(z_value), ncol = 5)
colnames(noise5) <- paste0("Noise", 1:5)
rownames(noise5) <- names(z_value)

################## chai #######################
set.seed(123)
chai_zhu_ae5 <- chai(z_value, ae5, B = 100)

length(lFDRselect(chai_zhu_ae5, 0.01, 1)) # 6
length(lFDRselect(chai_zhu_ae5, 0.05, 1)) # 38
length(lFDRselect(chai_zhu_ae5, 0.1, 1))  # 65


chai_zhu_pca5 <- chai(z_value, pca5, B = 100)

length(lFDRselect(chai_zhu_pca5, 0.01, 1)) # 6
length(lFDRselect(chai_zhu_pca5, 0.05, 1)) # 29
length(lFDRselect(chai_zhu_pca5, 0.1, 1))  # 55

chai_zhu_noise <- chai(z_value, noise, B = 100)
length(lFDRselect(chai_zhu_noise, 0.01, 1)) # 3
length(lFDRselect(chai_zhu_noise, 0.05, 1)) # 20
length(lFDRselect(chai_zhu_noise, 0.1, 1))  # 53

chai_zhu_noise2 <- chai(z_value, noise2, B = 100)

set.seed(123)
chai_zhu_noise5 <- chai(z_value, noise5, B = 100)
length(lFDRselect(chai_zhu_noise5, 0.01, 1))
length(lFDRselect(chai_zhu_noise5, 0.05, 1))
length(lFDRselect(chai_zhu_noise5, 0.1, 1))


length(intersect(which(p_adjusted < 0.05), lFDRselect(chai_zhu_ae5, 0.05, 1)))
length(intersect(which(p_adjusted < 0.1), lFDRselect(chai_zhu_ae5, 0.1, 1)))


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

chai_sel_ae5 <- label_direction(outcome, z_value, lFDRselect(chai_zhu_ae5, 0.05, 1))
# chai_sel$feature[chai_sel$direction == "higher_in_cirrhosis"]


# Visualization
chai_sel <- lFDRselect(chai_zhu_ae5, 0.05, 1)
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
ae5_sel <- lFDRselect(chai_zhu_ae5, 0.05, 1)
pca5_sel <- lFDRselect(chai_zhu_pca5, 0.05, 1)
noise1_sel <- lFDRselect(chai_zhu_noise, 0.05, 1)
noise2_sel <- lFDRselect(chai_zhu_noise2, 0.05, 1)
# Keep the Venn to 5 sets; include Efron in the UpSet below.

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
efron_sel_tax <- names(z_value)[efron]

all_species <- sort(unique(c(
  ae5_sel_tax,
  pca5_sel_tax,
  noise1_sel_tax,
  noise2_sel_tax,
  bh_sel_tax,
  efron_sel_tax
)))

upset_df <- tibble(species = all_species) %>%
  mutate(
    AE     = species %in% ae5_sel_tax,
    PCA    = species %in% pca5_sel_tax,
    noise1 = species %in% noise1_sel_tax,
    noise2 = species %in% noise2_sel_tax,
    BH     = species %in% bh_sel_tax,
    Efron_lFDR = species %in% efron_sel_tax
  )

upset_plot <- tryCatch(
  ComplexUpset::upset(
    upset_df,
    intersect = c("AE","PCA","noise1","noise2","BH","Efron_lFDR"),
    base_annotations = list(
      "Intersection size" = ComplexUpset::intersection_size())
  ),
  error = function(e) {
    message("Skipping UpSet plot: ", conditionMessage(e))
    NULL
  }
)

################## Adapt #######################
# alphas <- c(0.001, 0.01, 0.05)
alphas <- seq(0.01, 0.1, by = 0.01)


# GLM
# Builds one natural-spline formula per candidate df using the provided covariate names.
make_ns_formulas <- function(var_names, dfs) {
  lapply(dfs, function(d) {
    rhs <- paste0("ns(", var_names, ", df=", d, ")", collapse = " + ")
    as.formula(paste("~", rhs))
  })
}

# ns - AE5
x_ae5 <- data.frame(ae5)
dfs <- c(2, 4, 6, 8)
formula_ns_ae5 <- make_ns_formulas(colnames(x_ae5), dfs)

glm_ns_ae5 <- adapt_glm(x_ae5, pvals = p_value, alphas=alphas,
                        pi_formulas = formula_ns_ae5, mu_formulas = formula_ns_ae5)
# alpha = 0.1: FDPhat 0.0882, Number of Rej. 34
# alpha = 0.09: FDPhat 0.0882, Number of Rej. 34
# alpha = 0.08: FDPhat 0.069, Number of Rej. 29
# alpha = 0.07: FDPhat 0.069, Number of Rej. 29

length(intersect(which(p_adjusted <= 0.1), glm_ns_ae5$rejs[10][[1]]))

# ns - PCA5
x_pca5 <- data.frame(pca5)
dfs <- c(2, 4, 6, 8)
formula_ns_pca5 <- make_ns_formulas(colnames(x_pca5), dfs)

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

# ns - noise5
x_noise5 <- data.frame(noise5)
formula_ns_noise5 <- make_ns_formulas(colnames(x_noise5), dfs)

glm_ns_noise5 <- adapt_glm(x_noise5, pvals = p_value, alphas=alphas,
                           pi_formulas = formula_ns_noise5, mu_formulas = formula_ns_noise5)


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

# ns - noise5
gmm_ns_p_noise5 <- adapt_gmm(x_noise5, pvals = p_value,  alphas=alphas,
                             beta_formulas = formula_ns_noise5)

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

# ns - noise5
gmm_ns_z_noise5 <- adapt_gmm(x_noise5, z = z_value,  alphas=alphas,
                             beta_formulas = formula_ns_noise5, testing = "two_sided")

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
if (fdrreg_available) {
  set.seed(123)
  fdr_theo_ae5 <- FDRreg::FDRreg(z_value, as.matrix(ae5), nulltype = 'theoretical')
  length(which(fdr_theo_ae5$FDR <= 0.01))   # 2
  length(which(fdr_theo_ae5$FDR <= 0.05))   # 25
  length(which(fdr_theo_ae5$FDR <= 0.1))    # 58

  length(intersect(which(p_adjusted <= alphas[5]),  which(fdr_theo_ae5$FDR <= alphas[5])))

  fdr_theo_pca5 <- FDRreg::FDRreg(z_value, as.matrix(pca5), nulltype = 'theoretical')
  length(which(fdr_theo_pca5$FDR <= 0.01))   # 3
  length(which(fdr_theo_pca5$FDR <= 0.05))   # 23
  length(which(fdr_theo_pca5$FDR <= 0.1))    # 58

  fdr_theo_noise <- FDRreg::FDRreg(z_value, as.matrix(noise), nulltype = 'theoretical')
  length(which(fdr_theo_noise$FDR <= 0.01))   # 2
  length(which(fdr_theo_noise$FDR <= 0.05))   # 19
  length(which(fdr_theo_noise$FDR <= 0.1))    # 55

  fdr_theo_noise2 <- FDRreg::FDRreg(z_value, as.matrix(noise2), nulltype = 'theoretical')
  length(which(fdr_theo_noise2$FDR <= 0.01))   # 2
  length(which(fdr_theo_noise2$FDR <= 0.05))   # 18
  length(which(fdr_theo_noise2$FDR <= 0.1))    # 55

  fdr_theo_noise5 <- FDRreg::FDRreg(z_value, as.matrix(noise5), nulltype = 'theoretical')
  length(which(fdr_theo_noise5$FDR <= 0.01))
  length(which(fdr_theo_noise5$FDR <= 0.05))
  length(which(fdr_theo_noise5$FDR <= 0.1))
}

###########################################################
# Summarize into a different q level table
q_levels <- alphas

# chai
chai_rejs_ae5 <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_zhu_ae5, q, 1))
})

chai_rejs_pca5 <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_zhu_pca5, q, 1))
})

chai_rejs_noise <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_zhu_noise, q, 1))
})

chai_rejs_noise2 <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_zhu_noise2, q, 1))
})

chai_rejs_noise5 <- sapply(q_levels, function(q) {
  length(lFDRselect(chai_zhu_noise5, q, 1))
})

# adapt_glm
adapt_glm_ns_ae5 <- glm_ns_ae5$nrejs
adapt_glm_ns_pca5 <- glm_ns_pca5$nrejs
adapt_glm_ns_noise <- glm_ns_noise$nrejs
adapt_glm_ns_noise2 <- glm_ns_noise2$nrejs
adapt_glm_ns_noise5 <- glm_ns_noise5$nrejs


# adapt_gmm_p
adapt_gmm_p_ns_ae5 <- gmm_ns_p_ae5$nrejs
adapt_gmm_p_ns_pca5 <- gmm_ns_p_pca5$nrejs
adapt_gmm_p_ns_noise <- gmm_ns_p_noise$nrejs
adapt_gmm_p_ns_noise2 <- gmm_ns_p_noise2$nrejs
adapt_gmm_p_ns_noise5 <- gmm_ns_p_noise5$nrejs


# adapt_gmm_z
adapt_gmm_z_ns_ae5 <- gmm_ns_z_ae5$nrejs
adapt_gmm_z_ns_pca5 <- gmm_ns_z_pca5$nrejs
adapt_gmm_z_ns_noise <- gmm_ns_z_noise$nrejs
adapt_gmm_z_ns_noise2 <- gmm_ns_z_noise2$nrejs
adapt_gmm_z_ns_noise5 <- gmm_ns_z_noise5$nrejs


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
if (fdrreg_available) {
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

  rdrreg_noise5 <- sapply(q_levels, function(q) {
    length(which(fdr_theo_noise5$FDR <= q))
  })
}

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

# Two-stage Benjamini-Krieger-Yekutieli adaptive BH.
adaptive_bh_reject_count <- function(pvals, alpha) {
  m <- length(pvals)
  alpha_stage1 <- alpha / (1 + alpha)
  r1 <- sum(p.adjust(pvals, method = "BH") <= alpha_stage1)
  m0_hat <- max(1, m - r1)
  alpha_stage2 <- min(1, alpha * m / m0_hat)
  sum(p.adjust(pvals, method = "BH") <= alpha_stage2)
}

adaptive_bh_rejs <- sapply(q_levels, function(q) {
  adaptive_bh_reject_count(p_value, q)
})

# Storey q-value baseline estimated directly from the p-values.
estimate_storey_qvalues <- function(pvals, lambda = seq(0.05, 0.90, by = 0.05), spline_df = 3) {
  m <- length(pvals)
  pi0_lambda <- sapply(lambda, function(l) {
    mean(pvals > l) / (1 - l)
  })
  pi0_lambda <- pmin(pmax(pi0_lambda, 0), 1)
  spline_fit <- smooth.spline(lambda, pi0_lambda, df = spline_df)
  pi0 <- predict(spline_fit, x = max(lambda))$y
  pi0 <- min(max(pi0, 1 / m), 1)

  o <- order(pvals)
  p_sorted <- pvals[o]
  q_sorted <- numeric(m)
  q_sorted[m] <- pi0 * p_sorted[m]
  if (m > 1) {
    for (i in seq.int(m - 1, 1)) {
      q_sorted[i] <- min(pi0 * m * p_sorted[i] / i, q_sorted[i + 1])
    }
  }
  q_sorted <- pmin(pmax(q_sorted, 0), 1)

  qvals <- numeric(m)
  qvals[o] <- q_sorted
  qvals
}

storey_qvalues <- estimate_storey_qvalues(p_value)
storey_rejs <- sapply(q_levels, function(q) {
  sum(storey_qvalues <= q)
})

# Efron lFDR under theoretical N(0,1) null
efron_lfdr_rejs <- sapply(q_levels, function(q) {
  length(which(efron_lfdr <= q))
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
  IHW_x1 = ihw_x1_ae5,
  IHW_x2 = ihw_x2_ae5,
  IHW_x3 = ihw_x3_ae5,
  IHW_x4 = ihw_x4_ae5,
  IHW_x5 = ihw_x5_ae5,
  BH = bh_rejs,
  Efron_lFDR = efron_lfdr_rejs
)
if (fdrreg_available) {
  zhu_sum_ae5$FDRreg <- rdrreg_ae5
  zhu_sum_ae5 <- zhu_sum_ae5 %>% relocate(FDRreg, .after = OrderShapeEM_x5)
}

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
  IHW_x1 = ihw_x1_pca5,
  IHW_x2 = ihw_x2_pca5,
  IHW_x3 = ihw_x3_pca5,
  IHW_x4 = ihw_x4_pca5,
  IHW_x5 = ihw_x5_pca5,
  BH = bh_rejs,
  Efron_lFDR = efron_lfdr_rejs
)
if (fdrreg_available) {
  zhu_sum_pca5$FDRreg <- rdrreg_pca5
  zhu_sum_pca5 <- zhu_sum_pca5 %>% relocate(FDRreg, .after = OrderShapeEM_x5)
}


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
  mutate(Method = dplyr::recode(
    as.character(Method),
    Efron_lFDR = "Efron lFDR N(0,1)"
  ))

p_zhu_ae <- plot_rejections_vs_q(zhu_ae,
                              title = "Number of Discoveries vs q, with Autoencoder",
                              x_breaks = sort(unique(zhu_ae$q)))

pdf("./plots/Zhu/Zhu_AE_updated_locfdr.pdf", width = 10, height = 5)
plot(p_zhu_ae)
dev.off()

zhu_ae_updated <- zhu_ae %>%
  filter(Method != "Efron lFDR N(0,1)")

p_zhu_ae_updated <- plot_rejections_vs_q(
  zhu_ae_updated,
  title = "Number of Discoveries vs q, with Autoencoder",
  x_breaks = sort(unique(zhu_ae_updated$q)),
  label_nudge_multipliers = c("^OrderShapeEM" = 1.7)
)

pdf("./figures/Zhu_AE_updated_orderlabels_right.pdf", width = 10, height = 5)
plot(p_zhu_ae_updated)
dev.off()


# PCA
zhu_pca <- zhu_sum_pca5 %>%
  pivot_longer(
    cols = -q,
    names_to = "Method",
    values_to = "Rejections"
  ) %>%
  mutate(Method = dplyr::recode(
    as.character(Method),
    Efron_lFDR = "Efron lFDR N(0,1)"
  ))

p_zhu_pca <- plot_rejections_vs_q(zhu_pca,
                                 title = "Number of Discoveries vs q, with PCA",
                                 x_breaks = sort(unique(zhu_pca$q)))

pdf("./plots/Zhu/Zhu_PCA_updated_locfdr.pdf", width = 10, height = 5)
plot(p_zhu_pca)
dev.off()



######################################### AE vs. PCA #################################################
noise_set_labels <- sprintf("Noise%02d", seq_len(10))

generate_noise5_matrix <- function(seed) {
  set.seed(seed)
  noise_mat <- matrix(rnorm(length(z_value) * 5), nrow = length(z_value), ncol = 5)
  colnames(noise_mat) <- paste0("Noise", seq_len(5))
  rownames(noise_mat) <- names(z_value)
  noise_mat
}

fit_noise5_methods <- function(noise_mat) {
  x_noise <- data.frame(noise_mat)
  formula_ns_noise <- make_ns_formulas(colnames(x_noise), dfs)

  set.seed(123)
  chai_fit <- chai(z_value, noise_mat, B = 100)

  adapt_glm_fit <- adapt_glm(
    x_noise,
    pvals = p_value,
    alphas = alphas,
    pi_formulas = formula_ns_noise,
    mu_formulas = formula_ns_noise
  )

  adapt_gmm_p_fit <- adapt_gmm(
    x_noise,
    pvals = p_value,
    alphas = alphas,
    beta_formulas = formula_ns_noise
  )

  adapt_gmm_z_fit <- adapt_gmm(
    x_noise,
    z = z_value,
    alphas = alphas,
    beta_formulas = formula_ns_noise,
    testing = "two_sided"
  )

  result <- list(
    chai = chai_fit,
    adapt_glm = adapt_glm_fit,
    adapt_gmm_p = adapt_gmm_p_fit,
    adapt_gmm_z = adapt_gmm_z_fit
  )

  if (fdrreg_available) {
    set.seed(123)
    result$fdrreg <- FDRreg::FDRreg(z_value, as.matrix(noise_mat), nulltype = "theoretical")
  }

  result
}

noise5_fits <- vector("list", length(noise_set_labels))
names(noise5_fits) <- noise_set_labels
noise5_fits[[1]] <- list(
  chai = chai_zhu_noise5,
  adapt_glm = glm_ns_noise5,
  adapt_gmm_p = gmm_ns_p_noise5,
  adapt_gmm_z = gmm_ns_z_noise5
)
if (fdrreg_available) {
  noise5_fits[[1]]$fdrreg <- fdr_theo_noise5
}

for (i in 2:length(noise_set_labels)) {
  noise5_fits[[i]] <- fit_noise5_methods(generate_noise5_matrix(1000 + i))
}

zhu_sum_ae_pca <- data.frame(
  q = q_levels,
  chai_ae5 = chai_rejs_ae5,
  chai_pca5 = chai_rejs_pca5,
  adapt_glm_ae5 = adapt_glm_ns_ae5,
  adapt_glm_pca5 = adapt_glm_ns_pca5,
  adapt_gmm_p_ae5 = adapt_gmm_p_ns_ae5,
  adapt_gmm_p_pca5 = adapt_gmm_p_ns_pca5,
  adapt_gmm_z_ae5 = adapt_gmm_z_ns_ae5,
  adapt_gmm_z_pca5 = adapt_gmm_z_ns_pca5,
  BH = bh_rejs,
  Adaptive_BH = adaptive_bh_rejs,
  Storey = storey_rejs,
  Efron_lFDR = efron_lfdr_rejs
)

for (i in seq_along(noise_set_labels)) {
  noise_suffix <- tolower(noise_set_labels[[i]])
  fit_i <- noise5_fits[[i]]

  zhu_sum_ae_pca[[paste0("chai_", noise_suffix)]] <- sapply(q_levels, function(q) {
    length(lFDRselect(fit_i$chai, q, 1))
  })
  zhu_sum_ae_pca[[paste0("adapt_glm_", noise_suffix)]] <- fit_i$adapt_glm$nrejs
  zhu_sum_ae_pca[[paste0("adapt_gmm_p_", noise_suffix)]] <- fit_i$adapt_gmm_p$nrejs
  zhu_sum_ae_pca[[paste0("adapt_gmm_z_", noise_suffix)]] <- fit_i$adapt_gmm_z$nrejs
  if (fdrreg_available) {
    zhu_sum_ae_pca[[paste0("FDRreg_", noise_suffix)]] <- sapply(q_levels, function(q) {
      length(which(fit_i$fdrreg$FDR <= q))
    })
  }
}

if (fdrreg_available) {
  zhu_sum_ae_pca$FDRreg_ae5 <- rdrreg_ae5
  zhu_sum_ae_pca$FDRreg_pca5 <- rdrreg_pca5
}

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
      grepl("_ae5$", name) ~ "AE5",
      grepl("_pca5$", name) ~ "PCA5",
      grepl("_noise[0-9]{2}$", name) ~ sub(".*_(noise[0-9]{2})$", "\\1", name),
      TRUE                    ~ "baseline"
    ),
    Method = case_when(
      side == "baseline" ~ name,
      TRUE         ~ sub("_(ae5|pca5|noise[0-9]{2})$", "", name)
    )
  ) %>%
  mutate(
    side = ifelse(grepl("^noise", side), sub("^noise", "Noise", side), side)
  )

side_levels <- c("AE5", "PCA5", noise_set_labels)
zhu_long$side <- factor(zhu_long$side, levels = side_levels)

# Only methods with side information; baselines are added back into each facet.
zhu_long_methods <- zhu_long %>%
  filter(side != "baseline")

method_levels <- intersect(
  c("chai", "adapt_glm", "adapt_gmm_p", "adapt_gmm_z", "FDRreg"),
  unique(as.character(zhu_long_methods$Method))
)

zhu_long_methods <- zhu_long_methods %>%
  mutate(Method = factor(Method, levels = method_levels))

# --- End-of-line label positions (one per Method × side) -------------------
endpoints <- zhu_long_methods %>%
  arrange(Method, side, q) %>%
  group_by(Method, side) %>%
  slice_tail(n = 1) %>%
  ungroup()

side_label_endpoints <- endpoints %>%
  filter(side %in% c("AE5", "PCA5"))

dx <- diff(range(zhu_long_methods$q))

make_baseline_line <- function(method_name, label) {
  zhu_long %>%
    dplyr::filter(Method == method_name) %>%
    dplyr::select(q, Rejections) %>%
    mutate(Baseline = label) %>%
    tidyr::crossing(Method = method_levels) %>%
    mutate(Method = factor(Method, levels = method_levels))
}

baseline_lines <- bind_rows(
  make_baseline_line("BH", "BH"),
  make_baseline_line("Adaptive_BH", "Adaptive BH"),
  make_baseline_line("Storey", "Storey"),
  make_baseline_line("Efron_lFDR", "Efron lFDR N(0,1)")
)

baseline_endpoints <- baseline_lines %>%
  group_by(Method, Baseline) %>%
  slice_max(q, with_ties = FALSE) %>%
  ungroup()

p_ae_pca <- ggplot(zhu_long_methods,
                   aes(x = q, y = Rejections, color = side)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_text_repel(
    data = side_label_endpoints,
    aes(label = side),
    nudge_x      = dx * 0.03,
    hjust        = 0,
    direction    = "y",
    segment.size = 0.2,
    box.padding  = 0.1,
    point.padding= 0.1,
    show.legend  = FALSE
  ) +
  geom_line(
    data        = baseline_lines,
    aes(x = q, y = Rejections, group = Baseline, linetype = Baseline),
    inherit.aes = FALSE,
    color       = "black",
    linewidth   = 0.8
  ) +
  geom_text_repel(
    data        = baseline_endpoints,
    aes(x = q, y = Rejections, label = Baseline),
    inherit.aes = FALSE,
    nudge_x     = dx * 0.03,
    hjust       = 0,
    direction   = "y",
    segment.size = 0.2,
    box.padding = 0.1,
    point.padding = 0.1,
    size        = 3,
    show.legend = FALSE,
    color       = "black"
  ) +
  scale_linetype_manual(
    values = c(
      "BH" = "dashed",
      "Adaptive BH" = "longdash",
      "Storey" = "dotted",
      "Efron lFDR N(0,1)" = "dotdash"
    )
  ) +
  facet_wrap(~ Method) +  # <-- remove scales = "free_y"
  scale_x_continuous(breaks = q_levels) +
  scale_y_continuous(
    limits = c(
      0,
      max(c(zhu_long_methods$Rejections, baseline_lines$Rejections))
    )
  ) +
  scale_color_manual(
    values = c(
      AE5 = "red",
      PCA5 = "blue",
      setNames(gray.colors(length(noise_set_labels), start = 0.25, end = 0.85), noise_set_labels)
    )
  ) +
  labs(
    title  = "Number of Discoveries at different q level",
    x      = "q level",
    y      = "Number of Rejections",
    color  = "Side information: "
  ) +
  guides(
    color = guide_legend(ncol = 2),
    linetype = "none"
  ) +
  coord_cartesian(clip = "off") +
  expand_limits(x = max(zhu_long_methods$q) + dx * 0.05) +
  theme_minimal() +
  theme(
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin     = margin(5.5, 60, 5.5, 5.5)
  )


pdf("./plots/Zhu/Zhu_AE5_vs_PCA5_vs_noise_locfdr.pdf",
    width = 11, height = 7)
plot(p_ae_pca)
dev.off()
