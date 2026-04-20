# install.packages(c("readr","dplyr","stringr","tidyr"))  # if needed
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# Path to your file
path <- "results_melanoma_taxa_with_effects.csv"

# Load
df <- read_csv(path, show_col_types = FALSE)

# Treat every column except 'taxon' as result text
result_cols <- setdiff(names(df), "taxon")

pmid_tbl <- df %>%
  select(any_of(result_cols)) %>%
  pivot_longer(everything(), names_to = "field", values_to = "value") %>%
  filter(!is.na(value), value != "NA", nzchar(value)) %>%
  # Extract both PMID and APMID forms; capture the digits
  mutate(pmids = str_extract_all(value, "(?i)\\bA?PMID\\s*(\\d+)\\b")) %>%
  unnest(pmids, keep_empty = FALSE) %>%
  transmute(field, pmid = str_extract(pmids, "\\d+"))

# Counts
total_pmid_mentions <- nrow(pmid_tbl)
unique_pmids_count  <- n_distinct(pmid_tbl$pmid)

cat("Total PMID mentions across all cells: ", total_pmid_mentions, "\n", sep = "")
cat("Unique PMIDs: ", unique_pmids_count, "\n", sep = "")

# Optional: per-column breakdown
per_field <- pmid_tbl %>%
  group_by(field) %>%
  summarise(
    total_mentions = n(),
    unique_pmids   = n_distinct(pmid),
    .groups = "drop"
  )

print(per_field)

# Optional: save the unique PMIDs to a CSV for auditing
unique_pmids <- pmid_tbl %>% distinct(pmid) %>% arrange(pmid)
# write_csv(unique_pmids, "unique_pmids_from_results.csv")