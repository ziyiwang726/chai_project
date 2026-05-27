args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_file <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
if (is.na(script_file)) {
  script_file <- "scripts/export_zhu_chai_ae5_selected_species.R"
}
project_root <- normalizePath(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), ".."), mustWork = TRUE)
setwd(project_root)

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
  library(coin)
  library(mclust)
  library(admix)
  library(mvtnorm)
})

extract_species <- function(taxon) {
  sp <- ifelse(
    grepl("s__", taxon, fixed = TRUE),
    sub(".*s__([^|;]+).*", "\\1", taxon),
    taxon
  )
  sp <- gsub("^s__", "", sp)
  sp <- gsub("_", " ", sp)
  trimws(sp)
}

label_direction <- function(outcome, z_value, selected_idx) {
  first_level <- levels(outcome)[1]
  lab_hi <- paste0("higher_in_", first_level)
  lab_lo <- paste0("lower_in_", first_level)
  ifelse(z_value[selected_idx] > 0, lab_hi, lab_lo)
}

source_first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (!length(existing)) {
    stop("Could not locate required helper file. Tried: ", paste(paths, collapse = ", "))
  }
  source(existing[[1]])
}

source_first_existing(c("./code/chai.R", "./Gopalakrishnan/code/chai.R", "./Gopalakrishnan03252026/codeZW/chai.R"))
source_first_existing(c("./code/conditionalParam.R", "./Gopalakrishnan/code/conditionalParam.R", "./Gopalakrishnan03252026/codeZW/conditionalParam.R"))
source_first_existing(c("./code/naiveRemoveOneObs.R", "./Gopalakrishnan/code/naiveRemoveOneObs.R", "./Gopalakrishnan03252026/codeZW/naiveRemoveOneObs.R"))
source_first_existing(c("./code/rGaussianMix.R", "./Gopalakrishnan/code/rGaussianMix.R", "./Gopalakrishnan03252026/codeZW/rGaussianMix.R"))
source_first_existing(c("./Gopalakrishnan/code/utils.R", "./Gopalakrishnan03252026/codeZW/utils.R", "./code/utils.R"))

ae5 <- read.csv("./data/Zhu_autoencoder/outputs/ZhuF_2020_species_AE5.csv", row.names = 1, check.names = FALSE)

abun <- curatedMetagenomicData("ZhuF_2020.relative_abundance", dryrun = FALSE)
outcome <- as.factor(abun$`2021-03-31.ZhuF_2020.relative_abundance`[[5]])
abun_mat <- assay(abun[[1]], "relative_abundance")

taxa_names <- sub("^.*\\|", "", rownames(abun_mat))
rownames(abun_mat) <- extract_species(taxa_names)

common_species <- intersect(rownames(ae5), rownames(abun_mat))
abun_mat_common <- abun_mat[common_species, , drop = FALSE]

cs <- colSums(abun_mat_common, na.rm = TRUE)
cs[cs == 0] <- 1
abun_mat_norm <- sweep(abun_mat_common, 2, cs, FUN = "/")
t_abun <- t(abun_mat_norm)

p_z_values <- sapply(colnames(t_abun), function(feature) {
  features <- t_abun[, feature, drop = FALSE]
  wilcox <- coin::wilcox_test(features[, 1] ~ outcome)
  c(p_value = coin::pvalue(wilcox), z_value = coin::statistic(wilcox))
}, simplify = "matrix")
p_z_values <- t(p_z_values)

p_value <- p_z_values[, "p_value"]
z_value <- p_z_values[, "z_value"]
p_adjusted <- p.adjust(p_value, method = "BH")

set.seed(123)
chai_zhu_ae5 <- chai(z_value, ae5, B = 100)

selected_idx <- lFDRselect(chai_zhu_ae5, 0.05, 1)
selected_rank <- match(selected_idx, chai_zhu_ae5$ord)

selected <- data.frame(
  selection_rank = selected_rank,
  feature = names(z_value)[selected_idx],
  z_value = unname(z_value[selected_idx]),
  p_value = unname(p_value[selected_idx]),
  BH_adjusted_p_value = unname(p_adjusted[selected_idx]),
  chai_lFDR = unname(chai_zhu_ae5$clFDR[selected_idx]),
  chai_avgFDR = unname(chai_zhu_ae5$avgFDR[selected_rank]),
  direction = label_direction(outcome, z_value, selected_idx),
  ae5[selected_idx, , drop = FALSE],
  row.names = NULL,
  check.names = FALSE
)
selected <- selected[order(selected$selection_rank), , drop = FALSE]

out_file <- "./plots/Zhu/chai_selected_species_ae5_q005.csv"
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.csv(selected, out_file, row.names = FALSE)

message("Wrote ", nrow(selected), " selected species to ", normalizePath(out_file, mustWork = FALSE))
