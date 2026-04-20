## =========================================================
##  Setup and Load Data
## =========================================================
library(dplyr)
library(stringr)

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  idx <- grep("^--file=", cmd_args)
  if (length(idx)) return(normalizePath(sub("^--file=", "", cmd_args[idx[1]])))
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(normalizePath(ofile))
  NA_character_
}

script_path <- get_script_path()
script_dir <- if (!is.na(script_path)) dirname(script_path) else normalizePath(getwd())
llm_provider <- tolower(Sys.getenv("LLM_PROVIDER", "openai"))
provider_dir <- Sys.getenv(
  "SIDECOV_PROVIDER_DIR",
  switch(
    llm_provider,
    openai = "openAIGenerated",
    claude = "claudeGenerated",
    gemini = "geminiGenerated",
    combined = "combined",
    paste0(llm_provider, "Generated")
  )
)
llm_file_tag <- Sys.getenv("LLM_FILE_TAG", llm_provider)
sidecov_root <- file.path(script_dir, "..", "sideCov", provider_dir)

provider_tagged_file <- function(name) {
  if (!nzchar(llm_file_tag)) return(file.path(script_dir, name))
  sub("(\\.[^.]+)$", paste0("_", llm_file_tag, "\\1"), file.path(script_dir, name))
}

provider_tagged_basename <- function(name) {
  basename(provider_tagged_file(name))
}

score_input_file <- function(mode) {
  file.path(
    sidecov_root,
    mode,
    provider_tagged_basename(sprintf("taxon_ABCD_probit_stoufferZ_%s.csv", mode))
  )
}

cat("[CHECKPOINT] Loading taxonomy metadata files...\n")
tax_info <- read.csv(
  file.path(script_dir, "taxa_ids_filtered.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
d1_tax_raw <- read.csv(
  file.path(script_dir, "d1Taxonomy.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rank_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
rank_idx <- setNames(seq_along(rank_levels), rank_levels)
tax_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

clean_tax_string <- function(x) {
  x <- gsub("^[a-z]__", "", x, ignore.case = TRUE)
  x <- gsub("\\[|\\]", "", x)
  trimws(x)
}

d1_tax <- d1_tax_raw %>%
  mutate(across(all_of(tax_cols), clean_tax_string))

tax_info2 <- tax_info %>%
  mutate(
    taxon = clean_tax_string(taxon),
    rank_level = match(tolower(taxon_rank), rank_levels)
  )

score_cols <- c("A", "B", "C")
aux_rank_order <- c("phylum", "class", "order", "family", "genus", "species")

make_otu_table <- function(tsub, d1_tax, type = c("class", "order", "family", "genus", "species")) {
  type <- match.arg(type)

  lvl_seq <- switch(
    type,
    class = c("species", "genus", "family", "order", "class", "phylum"),
    order = c("species", "genus", "family", "order", "class", "phylum"),
    family = c("species", "genus", "family", "order", "class", "phylum"),
    genus = c("species", "genus", "family", "order", "class", "phylum"),
    species = c("species", "genus", "family", "order", "class", "phylum")
  )

  lookups <- lapply(score_cols, function(col_name) setNames(tsub[[col_name]], tsub$taxon))
  names(lookups) <- score_cols

  get_val <- function(row_idx, col_name) {
    for (lvl in lvl_seq) {
      tax_name <- d1_tax[[lvl]][row_idx]
      if (!is.na(tax_name) && nzchar(tax_name) && tax_name %in% names(lookups[[col_name]])) {
        val <- lookups[[col_name]][tax_name]
        if (!is.na(val)) return(val)
      }
    }
    0
  }

  n_rows <- nrow(d1_tax)
  out <- d1_tax %>% select(OTUname, all_of(tax_cols))
  for (col_name in score_cols) {
    out[[col_name]] <- vapply(seq_len(n_rows), get_val, numeric(1), col_name = col_name)
  }
  out
}

print_zero_counts <- function(df, name) {
  counts <- colSums(df[, score_cols, drop = FALSE] == 0, na.rm = TRUE)
  cat("\n--- Zero Counts for", name, "---\n")
  print(counts)
}

selection_target_folder <- function(target) {
  switch(
    tolower(target),
    family = "selectionTargetFamily",
    genus = "selectionTargetGenus",
    species = "selectionTargetSpecies",
    otu = "selectionTargetOTU",
    stop("Unsupported target: ", target)
  )
}

normalize_target_rank <- function(target) {
  if (tolower(target) == "otu") "species" else tolower(target)
}

is_valid_aux_combo <- function(target, aux_source) {
  target_idx <- match(normalize_target_rank(target), aux_rank_order)
  aux_idx <- match(tolower(aux_source), aux_rank_order)
  !is.na(target_idx) && !is.na(aux_idx) && aux_idx <= target_idx
}

search_ranks <- function(target, aux_source) {
  target_rank <- normalize_target_rank(target)
  target_idx <- match(target_rank, aux_rank_order)
  aux_idx <- match(tolower(aux_source), aux_rank_order)

  if (is.na(target_idx) || is.na(aux_idx)) return(character())

  # "aux_source and below" means the target rank plus only the allowed
  # coarser fallback ranks up to aux_source. If aux_source is finer than the
  # target rank, keep the search at the target rank only.
  if (aux_idx > target_idx) {
    return(target_rank)
  }

  rev(aux_rank_order[aux_idx:target_idx])
}

build_rank_lookups <- function(tab_ranked) {
  lookup_ranks <- c("phylum", "class", "order", "family", "genus", "species")
  out <- setNames(vector("list", length(score_cols)), score_cols)

  for (col_name in names(out)) {
    out[[col_name]] <- setNames(vector("list", length(lookup_ranks)), lookup_ranks)
    for (rank_name in lookup_ranks) {
      rank_tbl <- tab_ranked %>%
        filter(rank_level == rank_idx[rank_name]) %>%
        transmute(
          taxon = clean_tax_string(taxon),
          value = as.numeric(.data[[col_name]])
        ) %>%
        distinct(taxon, .keep_all = TRUE)
      out[[col_name]][[rank_name]] <- setNames(rank_tbl$value, rank_tbl$taxon)
    }
  }

  out
}

base_target_table <- function(target) {
  target_norm <- tolower(target)

  if (target_norm == "otu") {
    return(d1_tax %>% select(OTUname, all_of(tax_cols)))
  }

  keep_cols <- tax_cols[seq_len(match(target_norm, tax_cols))]
  d1_tax %>%
    select(all_of(keep_cols)) %>%
    filter(!is.na(.data[[target_norm]]), .data[[target_norm]] != "") %>%
    distinct() %>%
    arrange(.data[[target_norm]])
}

fallback_value <- function(row_df, col_name, ranks_to_search, lookups) {
  for (rank_name in ranks_to_search) {
    tax_name <- row_df[[rank_name]][[1]]
    if (!is.na(tax_name) && nzchar(tax_name)) {
      val <- lookups[[col_name]][[rank_name]][tax_name]
      if (length(val) == 1 && !is.na(val)) {
        return(as.numeric(val))
      }
    }
  }
  0
}

make_target_table <- function(tab_ranked, target, aux_source) {
  ranks_to_search <- search_ranks(target, aux_source)
  if (!length(ranks_to_search)) return(NULL)

  out <- base_target_table(target)
  lookups <- build_rank_lookups(tab_ranked)

  for (col_name in score_cols) {
    out[[col_name]] <- vapply(
      seq_len(nrow(out)),
      function(i) fallback_value(out[i, , drop = FALSE], col_name, ranks_to_search, lookups),
      numeric(1)
    )
  }

  out
}

write_sidecov_table <- function(df, suffix, target, aux_source) {
  target_folder <- selection_target_folder(target)
  out_dir <- file.path(sidecov_root, suffix, target_folder)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  target_label <- if (tolower(target) == "otu") "OTU" else tolower(target)
  provider_suffix <- if (identical(llm_provider, "openai") && nzchar(llm_file_tag)) {
    paste0("_", llm_file_tag)
  } else {
    ""
  }
  out_file <- file.path(
    out_dir,
    paste0("d1Taxonomy_", target_label, "_ABCD_", aux_source, "_and_below_", suffix, provider_suffix, ".csv")
  )

  write.csv(df, out_file, row.names = FALSE)
  cat("  -", out_file, "\n")
}

write_all_sidecov_tables <- function(tab_ranked, suffix) {
  targets <- c("family", "genus", "species", "otu")
  aux_sources <- c("phylum", "class", "order", "family", "genus", "species")

  cat("\n[CHECKPOINT] Writing sideCov outputs for", suffix, "to", sidecov_root, "...\n")

  for (target in targets) {
    for (aux_source in aux_sources) {
      if (!is_valid_aux_combo(target, aux_source)) next
      tbl <- make_target_table(tab_ranked, target, aux_source)
      if (is.null(tbl)) next
      write_sidecov_table(tbl, suffix, target, aux_source)
    }
  }
}

run_mode <- function(input_csv, suffix) {
  if (!file.exists(input_csv)) {
    cat("[CHECKPOINT] Missing input, skipping:", input_csv, "\n")
    return(invisible(NULL))
  }

  cat("\n[CHECKPOINT] Processing mode:", suffix, "using", input_csv, "\n")
  tab_abcd <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)

  tab_ranked <- tab_abcd %>%
    left_join(tax_info2 %>% select(taxon, rank_level), by = "taxon")

  tab_class <- tab_ranked %>%
    filter(!is.na(rank_level), rank_level >= rank_idx["phylum"]) %>%
    select(taxon, all_of(score_cols))

  tab_order <- tab_ranked %>%
    filter(!is.na(rank_level), rank_level >= rank_idx["order"]) %>%
    select(taxon, all_of(score_cols))

  tab_family <- tab_ranked %>%
    filter(!is.na(rank_level), rank_level >= rank_idx["family"]) %>%
    select(taxon, all_of(score_cols))

  tab_genus <- tab_ranked %>%
    filter(!is.na(rank_level), rank_level >= rank_idx["genus"]) %>%
    select(taxon, all_of(score_cols))

  tab_species <- tab_ranked %>%
    filter(!is.na(rank_level), rank_level >= rank_idx["species"]) %>%
    select(taxon, all_of(score_cols))

  otu_class <- make_otu_table(tab_class, d1_tax, type = "class")
  otu_order <- make_otu_table(tab_order, d1_tax, type = "order")
  otu_family <- make_otu_table(tab_family, d1_tax, type = "family")
  otu_genus <- make_otu_table(tab_genus, d1_tax, type = "genus")
  otu_species <- make_otu_table(tab_species, d1_tax, type = "species")

  print_zero_counts(otu_class, paste("Class & below", suffix))
  print_zero_counts(otu_order, paste("Order & below", suffix))
  print_zero_counts(otu_family, paste("Family & below", suffix))
  print_zero_counts(otu_genus, paste("Genus & below", suffix))
  print_zero_counts(otu_species, paste("Species & below", suffix))
  cat("[CHECKPOINT] OTU-level outputs will be written only under sideCov selectionTargetOTU.\n")

  write_all_sidecov_tables(tab_ranked, suffix)
}

run_mode(score_input_file("unweighted"), "unweighted")
run_mode(score_input_file("weighted"), "weighted")

cat("\n[CHECKPOINT] produceFinalTable complete.\n")
