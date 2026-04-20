# Run this script ONCE in a fresh R session before running real_data_16s_LLM_plot.R
# In RStudio: Session -> Restart R, then source this file.

# macOS toolchain needed by some source packages (e.g., FDRreg)
# On macOS, install the GNU Fortran compiler first: https://mac.r-project.org/tools/
flibs_from_r <- function() {
  r_bin <- file.path(R.home("bin"), "R")
  out <- tryCatch(
    system2(r_bin, c("CMD", "config", "FLIBS"), stdout = TRUE, stderr = FALSE),
    warning = function(w) character(),
    error = function(e) character()
  )
  paste(out, collapse = " ")
}

fortran_runtime_ready <- function() {
  flibs <- flibs_from_r()
  if (!nzchar(flibs))
    return(nzchar(Sys.which("gfortran")))

  lib_dirs <- unique(sub("^-L", "", unlist(regmatches(
    flibs, gregexpr("-L[^[:space:]]+", flibs)
  ))))
  dirs_ok <- length(lib_dirs) > 0 && all(dir.exists(lib_dirs))
  has_emutls <- !grepl("-lemutls_w", flibs, fixed = TRUE) ||
    any(file.exists(file.path(lib_dirs, c("libemutls_w.a", "libemutls_w.dylib"))))
  has_heapt <- !grepl("-lheapt_w", flibs, fixed = TRUE) ||
    any(file.exists(file.path(lib_dirs, c("libheapt_w.a", "libheapt_w.dylib"))))
  dirs_ok && has_emutls && has_heapt
}

install_macos_gfortran <- function() {
  result <- list(ready = TRUE, manual_cmd = NA_character_)
  if (Sys.info()[["sysname"]] != "Darwin" || fortran_runtime_ready())
    return(result)

  tools_url <- "https://mac.r-project.org/tools/"
  page <- tryCatch(system2("curl", c("-fsSL", tools_url), stdout = TRUE),
                   warning = function(w) character(),
                   error = function(e) character())
  pkg_candidates <- unique(unlist(regmatches(
    page, gregexpr("gfortran-[0-9][^\"'[:space:]]*\\.pkg", page)
  )))

  if (!length(pkg_candidates)) {
    result$ready <- FALSE
    message("Could not find a gfortran installer on ", tools_url,
            ". Please install it manually, then rerun this script.")
    return(result)
  }

  versions <- sub("^gfortran-([0-9]+(\\.[0-9]+)*).*$", "\\1", pkg_candidates)
  latest_pkg <- pkg_candidates[order(numeric_version(versions), decreasing = TRUE)][1]
  pkg_url <- paste0(tools_url, latest_pkg)
  pkg_path <- file.path(tempdir(), latest_pkg)
  manual_path <- file.path("/tmp", latest_pkg)
  result$manual_cmd <- paste(
    "curl -fL", shQuote(pkg_url), "-o", shQuote(manual_path),
    "&& sudo installer -pkg", shQuote(manual_path), "-target /"
  )

  auto_install <- isTRUE(getOption("install_packages.auto_install_gfortran", FALSE))
  if (!auto_install) {
    result$ready <- FALSE
    message("GNU Fortran runtime is missing.")
    message("Run this command in Terminal, then rerun this script:\n", result$manual_cmd)
    message("Tip: set options(install_packages.auto_install_gfortran = TRUE) ",
            "before sourcing to let this script try auto-install.")
    return(result)
  }

  if (system2("curl", c("-fL", "-o", pkg_path, pkg_url)) != 0L) {
    result$ready <- FALSE
    message("Failed to download gfortran from ", pkg_url,
            ". Run this command manually:\n", result$manual_cmd)
    return(result)
  }

  if (system2("sudo", c("installer", "-pkg", pkg_path, "-target", "/")) != 0L) {
    result$ready <- FALSE
    message("gfortran installation failed. Run this command manually:\n", result$manual_cmd)
    return(result)
  }

  result$ready <- fortran_runtime_ready()
  if (!result$ready)
    message("GNU Fortran installation appears incomplete. Run this command manually:\n", result$manual_cmd)
  result
}

gfortran_setup <- install_macos_gfortran()
failed_pkgs <- character()

# CRAN packages
cran_pkgs <- c(
  "devtools", "ape", "coin", "locfdr", "admix", "mvtnorm", "mixtools",
  "KScorrect", "mclust", "splines2", "adaptMT", "Iso",
  "dplyr", "tibble", "stringr", "tidyr", "readr", "auctestr",
  "GGally", "ggplot2", "ggrepel", "viridisLite", "pheatmap", "ggnewscale"
)
to_install <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
bioc_pkgs <- c("phyloseq", "IHW", "ggtree", "ggtreeExtra", "qvalue")
to_install_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install_bioc)) BiocManager::install(to_install_bioc)

# GitHub packages
if (!requireNamespace("CATMicrobiome", quietly = TRUE))
  devtools::install_github("YushuShi/CATMicrobiome")
if (!requireNamespace("AdaPTGMM", quietly = TRUE))
  devtools::install_github("patrickrchao/AdaPTGMM")
if (!requireNamespace("adaptMT", quietly = TRUE))
  devtools::install_github("lihualei71/adaptMT")
  if (!requireNamespace("FDRreg", quietly = TRUE)) {
  if (Sys.info()[["sysname"]] == "Darwin" && !isTRUE(gfortran_setup$ready)) {
    failed_pkgs <- c(failed_pkgs, "FDRreg")
    message("Skipping FDRreg because GNU Fortran runtime is missing.")
    if (!is.na(gfortran_setup$manual_cmd) && nzchar(gfortran_setup$manual_cmd))
      message("Install it using:\n", gfortran_setup$manual_cmd)
  } else {
    tryCatch(
      devtools::install_github("cran/FDRreg"),
      error = function(e) {
        failed_pkgs <<- c(failed_pkgs, "FDRreg")
        message("FDRreg failed to install: ", conditionMessage(e))
        if (Sys.info()[["sysname"]] == "Darwin" &&
            !is.na(gfortran_setup$manual_cmd) && nzchar(gfortran_setup$manual_cmd))
          message("Install GNU Fortran from https://mac.r-project.org/tools/ and retry:\n",
                  gfortran_setup$manual_cmd)
      }
    )
  }
}
if (!requireNamespace("OrderShapeEM", quietly = TRUE))
  devtools::install_github("jchen1981/OrderShapeEM")

failed_pkgs <- unique(failed_pkgs)
if (length(failed_pkgs)) {
  message("Package installation completed with issues: ",
          paste(failed_pkgs, collapse = ", "),
          ". Install prerequisites and rerun this script.")
} else {
  message("All packages installed. Restart R before running real_data_16s_LLM_plot.R")
}
