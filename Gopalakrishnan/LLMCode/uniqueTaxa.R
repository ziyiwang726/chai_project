library(stringr)

args <- commandArgs(trailingOnly = TRUE)
in_csv <- if (length(args) >= 1) args[[1]] else "deepakTaxonomy.csv"
out_csv <- if (length(args) >= 2) args[[2]] else "deepakTaxaUnique.csv"

cat("[CHECKPOINT] Reading taxonomy table:", in_csv, "\n")
data <- read.csv(in_csv, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
data <- data[, -1, drop = FALSE]
dataVec <- as.vector(data)
dataVec <- unlist(dataVec)
dataVec <- unique(dataVec)
dataVec <- substring(dataVec, 4)
# dataVec: character vector
cleanVec <- dataVec |>
  str_remove_all("\\[|\\]") |>                                   # remove [ ]
  str_replace_all("_", " ") |>                                   # _ -> space
  str_remove_all("['\u2019]") |>                                 # remove ' and ’
  str_remove_all(regex("\\bunclassified\\b", ignore_case = TRUE))|> # drop "unclassified"
  str_squish()                                                   # trim & collapse spaces

# Split entries that have EXACTLY two words into two entries; keep others unchanged
outVec <- unlist(lapply(cleanVec, function(s) {
  if (!nzchar(s)) return(character(0))
  if (str_count(s, boundary("word")) == 2) {
    str_split(s, "\\s+", n = 2)[[1]]
  } else {
    s
  }
}), use.names = FALSE)                                # trim ends

outVec <- unique(outVec)

write.csv(outVec, file = out_csv, row.names = FALSE)
cat("[CHECKPOINT] Wrote unique taxa list:", out_csv, "\n")
