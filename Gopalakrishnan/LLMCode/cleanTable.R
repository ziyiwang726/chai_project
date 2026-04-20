## ---------------------------------------------------------------
##  Load packages
## ---------------------------------------------------------------
library(dplyr)
library(tidyr)
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
summary_dir <- file.path(sidecov_root, "summaryTables")

provider_tagged_file <- function(name) {
  if (!nzchar(llm_file_tag)) return(file.path(script_dir, name))
  sub("(\\.[^.]+)$", paste0("_", llm_file_tag, "\\1"), file.path(script_dir, name))
}

provider_tagged_basename <- function(name) {
  basename(provider_tagged_file(name))
}

clean_tax_string <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("^[a-z]__", "", x, ignore.case = TRUE)
  x <- gsub("\\[|\\]", "", x)
  x <- gsub("\\s+", " ", x)
  trimws(tolower(x))
}

escape_regex <- function(x) {
  str_replace_all(x, "([.|()\\^{}+$*?]|\\[|\\]|\\\\)", "\\\\\\1")
}

contains_name <- function(text, name) {
  if (is.na(text) || !nzchar(text) || is.na(name) || !nzchar(name)) {
    return(FALSE)
  }
  str_detect(
    text,
    regex(
      paste0("(^|[^[:alnum:]])", escape_regex(name), "($|[^[:alnum:]])"),
      ignore_case = TRUE
    )
  )
}

rank_levels <- c(
  "superkingdom", "domain", "kingdom", "phylum",
  "class", "order", "family", "genus", "species", "strain"
)
rank_idx <- setNames(seq_along(rank_levels), rank_levels)
lower_marker_regex <- regex(
  "(^|\\s)(otu|otus|asv|asvs|sgb|sgbs|strain|strains|sequence variant|sequence variants|amplicon|amplicons)(\\s|$)",
  ignore_case = TRUE
)
mouse_marker_regex <- regex("(^|\\s)(mouse|mice|murine)(\\s|$)", ignore_case = TRUE)

extract_species_genera <- function(text) {
  if (is.na(text) || !nzchar(text)) {
    return(character())
  }
  matches <- str_match_all(text, "\\b([A-Z][a-z]+)\\s+([a-z][A-Za-z0-9._-]+)\\b")[[1]]
  if (nrow(matches) == 0) {
    return(character())
  }
  genera <- tolower(matches[, 2])
  genera <- genera[!(tolower(matches[, 3]) %in% c("sp", "spp", "cf", "aff"))]
  unique(genera[nzchar(genera)])
}

tax_info <- read.csv(
  file.path(script_dir, "taxa_ids_filtered.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
tax_info_lookup <- tax_info %>%
  transmute(
    taxon_key = clean_tax_string(taxon),
    current_scientific_name = clean_tax_string(current_scientific_name),
    taxon_rank = clean_tax_string(taxon_rank)
  ) %>%
  distinct(taxon_key, .keep_all = TRUE)

known_taxa <- bind_rows(
  tax_info %>%
    transmute(name = clean_tax_string(taxon), rank = clean_tax_string(taxon_rank)),
  tax_info %>%
    transmute(name = clean_tax_string(current_scientific_name), rank = clean_tax_string(taxon_rank))
) %>%
  filter(name != "", !is.na(rank), rank %in% names(rank_idx)) %>%
  distinct()

cell_specificity_ok <- function(raw, taxon_name, current_name, taxon_rank) {
  effect_text <- str_match(raw, "\\[(.*)\\]\\s*$")[, 2]
  effect_norm <- clean_tax_string(effect_text)
  if (!nzchar(effect_norm)) {
    return(TRUE)
  }
  if (str_detect(effect_norm, mouse_marker_regex)) {
    return(FALSE)
  }

  allowed_names <- unique(Filter(nzchar, c(
    clean_tax_string(taxon_name),
    clean_tax_string(current_name)
  )))
  rank_norm <- clean_tax_string(taxon_rank)
  species_genera <- extract_species_genera(effect_text)

  if (identical(rank_norm, "species")) {
    other_species <- known_taxa$name[
      known_taxa$rank == "species" &
      !(known_taxa$name %in% allowed_names)
    ]
    if (any(vapply(other_species, function(name) contains_name(effect_norm, name), logical(1)))) {
      return(FALSE)
    }
    if (str_detect(effect_norm, lower_marker_regex)) {
      return(any(vapply(allowed_names, function(name) contains_name(effect_norm, name), logical(1))))
    }
    return(TRUE)
  }

  if (identical(rank_norm, "genus")) {
    if (length(species_genera)) {
      return(FALSE)
    }
    if (str_detect(effect_norm, lower_marker_regex)) {
      return(any(vapply(allowed_names, function(name) contains_name(effect_norm, name), logical(1))))
    }
  }

  if (length(species_genera)) {
    return(FALSE)
  }
  if (str_detect(effect_norm, lower_marker_regex)) {
    return(FALSE)
  }

  target_idx <- unname(rank_idx[[rank_norm]])
  if (is.null(target_idx) || is.na(target_idx)) {
    return(TRUE)
  }

  finer_names <- known_taxa$name[
    !(known_taxa$name %in% allowed_names) &
    !is.na(rank_idx[known_taxa$rank]) &
    rank_idx[known_taxa$rank] > target_idx
  ]
  !any(vapply(finer_names, function(name) contains_name(effect_norm, name), logical(1)))
}

input_matrix <- Sys.getenv("LLM_FILLED_FIXED_FILE", provider_tagged_file("taxon_article_matrix_filled.csv"))

cat("[CHECKPOINT] Loading", basename(input_matrix), "\n")
dat <- read.csv(
  input_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## ---------------------------------------------------------------
##  Reshape to long format
## ---------------------------------------------------------------
long <- dat %>%
  pivot_longer(
    cols = -taxon,
    names_to = "article",
    values_to = "raw"
  ) %>%
  filter(!is.na(raw), raw != "")

cat("[CHECKPOINT] Non-empty taxon-article cells:", nrow(long), "\n")

## ---------------------------------------------------------------
##  Parse extraction cell
##  Expected form: +0.05 A {n=120} [effect text]
##  Rules:
##  - direction must be + or -
##  - if p missing from parser, default to 0.05
##  - unknown direction (?) is excluded
##  - legacy D is merged into C
## ---------------------------------------------------------------
parse_cell <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) {
    return(data.frame(
      sign_chr = NA_character_,
      p = NA_real_,
      group = NA_character_,
      n_obs = NA_real_,
      effect_text = NA_character_
    ))
  }

  s <- trimws(x)
  effect_text <- str_match(s, "\\[(.*)\\]\\s*$")[, 2]

  n1 <- str_match(s, "(?i)\\{\\s*n\\s*=\\s*([0-9]+)\\s*\\}")[, 2]
  n2 <- str_match(s, "(?i)\\bn\\s*[=:]\\s*([0-9]+)")[, 2]
  n_obs <- suppressWarnings(as.numeric(ifelse(!is.na(n1), n1, n2)))

  core <- gsub("\\s*\\[.*$", "", s)
  m <- str_match(core, "^\\s*([+-])\\s*([0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\\s*([ABCD])\\b")

  sign_chr <- m[, 2]
  p_val <- suppressWarnings(as.numeric(m[, 3]))
  grp <- m[, 4]
  if (!is.na(grp) && grp == "D") {
    grp <- "C"
  }

  if (is.na(p_val)) {
    m2 <- str_match(core, "^\\s*([+-])\\s*([ABCD])\\b")
    if (!is.na(m2[, 2]) && !is.na(m2[, 3])) {
      sign_chr <- m2[, 2]
      grp <- m2[, 3]
      if (!is.na(grp) && grp == "D") {
        grp <- "C"
      }
      p_val <- 0.05
    }
  }

  if (is.na(sign_chr) || is.na(grp)) {
    return(data.frame(
      sign_chr = NA_character_,
      p = NA_real_,
      group = NA_character_,
      n_obs = n_obs,
      effect_text = effect_text
    ))
  }

  data.frame(
    sign_chr = sign_chr,
    p = p_val,
    group = grp,
    n_obs = n_obs,
    effect_text = effect_text
  )
}

parsed <- bind_rows(lapply(long$raw, parse_cell))
long2 <- bind_cols(long %>% select(taxon, article, raw), parsed) %>%
  filter(!is.na(sign_chr), !is.na(group), !is.na(p)) %>%
  mutate(
    p = ifelse(is.finite(p) & p > 0 & p <= 1, p, 0.05),
    sign_num = ifelse(sign_chr == "+", 1, -1),
    weight = ifelse(!is.na(n_obs) & n_obs > 0, n_obs, 1),
    taxon_key = clean_tax_string(taxon)
  ) %>%
  left_join(tax_info_lookup, by = "taxon_key")

specificity_ok <- mapply(
  cell_specificity_ok,
  long2$raw,
  long2$taxon,
  long2$current_scientific_name,
  long2$taxon_rank,
  USE.NAMES = FALSE
)
removed_n <- sum(!specificity_ok)
long2 <- long2[specificity_ok, , drop = FALSE]

cat("[CHECKPOINT] Parsed valid directional significant cells:", nrow(long2), "\n")
cat("[CHECKPOINT] Removed cells failing taxonomic specificity:", removed_n, "\n")

## ---------------------------------------------------------------
##  Table 1: taxon + A/B/C signed p-values (comma-separated)
## ---------------------------------------------------------------
table1 <- long2 %>%
  mutate(
    value = sign_num * p,
    value_str = format(value, scientific = FALSE, trim = TRUE)
  ) %>%
  group_by(taxon, group) %>%
  summarise(values = paste(value_str, collapse = ", "), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = values)

for (cc in c("A", "B", "C")) {
  if (!cc %in% names(table1)) {
    table1[[cc]] <- NA_character_
  }
}

table1 <- table1 %>% select(taxon, A, B, C) %>% arrange(taxon)

## ---------------------------------------------------------------
##  Probit transform and Stouffer combination
##  z_i = sign * qnorm(1 - p_i/2)  [two-sided p to signed z]
## ---------------------------------------------------------------
eps <- 1e-15

long_z <- long2 %>%
  mutate(
    p_clamped = pmin(pmax(p, eps), 1 - eps),
    z = sign_num * qnorm(1 - p_clamped / 2)
  )

table_unweighted <- long_z %>%
  group_by(taxon, group) %>%
  summarise(
    z_stouffer = {
      k <- sum(!is.na(z))
      if (k > 0) sum(z, na.rm = TRUE) / sqrt(k) else NA_real_
    },
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = group, values_from = z_stouffer)

for (cc in c("A", "B", "C")) {
  if (!cc %in% names(table_unweighted)) {
    table_unweighted[[cc]] <- NA_real_
  }
}

table_unweighted <- table_unweighted %>% select(taxon, A, B, C) %>% arrange(taxon)

table_weighted <- long_z %>%
  group_by(taxon, group) %>%
  summarise(
    z_stouffer_weighted = {
      ok <- !is.na(z) & !is.na(weight) & is.finite(weight) & weight > 0
      if (!any(ok)) {
        NA_real_
      } else {
        wz <- sum(weight[ok] * z[ok])
        denom <- sqrt(sum(weight[ok]^2))
        if (denom > 0) wz / denom else NA_real_
      }
    },
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = group, values_from = z_stouffer_weighted)

for (cc in c("A", "B", "C")) {
  if (!cc %in% names(table_weighted)) {
    table_weighted[[cc]] <- NA_real_
  }
}

table_weighted <- table_weighted %>% select(taxon, A, B, C) %>% arrange(taxon)

## ---------------------------------------------------------------
##  Save outputs
## ---------------------------------------------------------------
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sidecov_root, "unweighted"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(sidecov_root, "weighted"), recursive = TRUE, showWarnings = FALSE)

out_table1 <- file.path(summary_dir, provider_tagged_basename("taxon_ABCD.csv"))
out_unweighted <- file.path(
  sidecov_root,
  "unweighted",
  provider_tagged_basename("taxon_ABCD_probit_stoufferZ_unweighted.csv")
)
out_weighted <- file.path(
  sidecov_root,
  "weighted",
  provider_tagged_basename("taxon_ABCD_probit_stoufferZ_weighted.csv")
)

write.csv(table1, out_table1, row.names = FALSE)
write.csv(table_unweighted, out_unweighted, row.names = FALSE)
write.csv(table_weighted, out_weighted, row.names = FALSE)

cat("[CHECKPOINT] Saved:\n")
cat("  -", out_table1, "\n")
cat("  -", out_unweighted, "\n")
cat("  -", out_weighted, "\n")
