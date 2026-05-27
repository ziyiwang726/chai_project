#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript export_feature_selection_combined_table.R <input_rdata> <sidecov_csv> <output_csv>")
}

input_rdata <- args[[1]]
sidecov_csv <- args[[2]]
output_csv <- args[[3]]

clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    str_replace("^.*?__", "") %>%
    str_replace_all("\\[|\\]", "") %>%
    str_replace_all("^'+|'+$", "") %>%
    str_squish()
}

load(input_rdata)

required_objects <- c(
  "selected_features_long",
  "chai_selected_q05",
  "chai_selected_q10",
  "bh_selected_q05",
  "bh_selected_q10"
)
missing_objects <- setdiff(required_objects, ls())
if (length(missing_objects)) {
  stop("Missing objects in RData: ", paste(missing_objects, collapse = ", "))
}

if (!exists("selection_metadata") || is.null(selection_metadata$target_rank)) {
  stop("RData does not contain selection_metadata$target_rank.")
}

target_rank <- tolower(selection_metadata$target_rank[[1]])

base_info <- selected_features_long %>%
  select(otu_rep, rank_value, taxon_key, p_main, z_main, xA, xB, xC, p_adjusted, clFDR) %>%
  distinct() %>%
  mutate(taxon_key_join = clean_taxon(taxon_key))

add_flag <- function(df, flag_name) {
  if (!nrow(df)) {
    return(tibble(taxon_key_join = character(0), !!flag_name := logical(0)))
  }

  df %>%
    transmute(taxon_key_join = clean_taxon(taxon_key), !!flag_name := TRUE) %>%
    distinct()
}

flag_table <- list(
  add_flag(chai_selected_q05, "chai_q05"),
  add_flag(chai_selected_q10, "chai_q10"),
  add_flag(bh_selected_q05, "bh_q05"),
  add_flag(bh_selected_q10, "bh_q10")
) %>%
  bind_rows() %>%
  group_by(taxon_key_join) %>%
  summarise(
    chai_q05 = any(ifelse(is.na(chai_q05), FALSE, chai_q05)),
    chai_q10 = any(ifelse(is.na(chai_q10), FALSE, chai_q10)),
    bh_q05 = any(ifelse(is.na(bh_q05), FALSE, bh_q05)),
    bh_q10 = any(ifelse(is.na(bh_q10), FALSE, bh_q10)),
    .groups = "drop"
  )

side_raw <- read_csv(sidecov_csv, show_col_types = FALSE)

target_rank_col <- if (target_rank == "otu") {
  otu_candidates <- names(side_raw)[tolower(names(side_raw)) %in% c("otuname", "otu", "otu_rep")]
  if (!length(otu_candidates)) {
    stop("Could not find an OTU identifier column in side-information CSV.")
  }
  otu_candidates[[1]]
} else {
  rank_matches <- names(side_raw)[tolower(names(side_raw)) == target_rank]
  if (!length(rank_matches)) {
    stop("Could not find target rank column '", target_rank, "' in side-information CSV.")
  }
  rank_matches[[1]]
}

taxonomy_cols <- intersect(
  c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
  names(side_raw)
)

side_info <- side_raw %>%
  mutate(taxon_key_join = clean_taxon(.data[[target_rank_col]])) %>%
  group_by(taxon_key_join) %>%
  summarise(
    matched_rows_in_sidecov = n(),
    across(
      any_of(taxonomy_cols),
      ~ paste(unique(.x[!is.na(.x) & .x != ""]), collapse = " | ")
    ),
    target_value_from_sidecov = paste(
      unique(.data[[target_rank_col]][!is.na(.data[[target_rank_col]]) & .data[[target_rank_col]] != ""]),
      collapse = " | "
    ),
    A_sidecov = paste(unique(format(signif(A, 8), scientific = FALSE, trim = TRUE)), collapse = " | "),
    B_sidecov = paste(unique(format(signif(B, 8), scientific = FALSE, trim = TRUE)), collapse = " | "),
    C_sidecov = paste(unique(format(signif(C, 8), scientific = FALSE, trim = TRUE)), collapse = " | "),
    .groups = "drop"
  )

out <- base_info %>%
  left_join(flag_table, by = "taxon_key_join") %>%
  left_join(side_info, by = "taxon_key_join") %>%
  mutate(
    chai_q05 = ifelse(is.na(chai_q05), FALSE, chai_q05),
    chai_q10 = ifelse(is.na(chai_q10), FALSE, chai_q10),
    bh_q05 = ifelse(is.na(bh_q05), FALSE, bh_q05),
    bh_q10 = ifelse(is.na(bh_q10), FALSE, bh_q10)
  ) %>%
  select(
    otu_rep,
    rank_value,
    taxon_key,
    any_of(taxonomy_cols),
    target_value_from_sidecov,
    matched_rows_in_sidecov,
    p_main,
    z_main,
    p_adjusted,
    clFDR,
    xA,
    xB,
    xC,
    A_sidecov,
    B_sidecov,
    C_sidecov,
    chai_q05,
    chai_q10,
    bh_q05,
    bh_q10
  ) %>%
  arrange(desc(chai_q10), desc(chai_q05), clFDR, p_main)

write_csv(out, output_csv)
message("Wrote combined table: ", output_csv)
