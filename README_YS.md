# chai Project Repository
**chai** - Conditional Hypothesis testing using Auxiliary Information

**chai** is a covariate-informed statistical framework. It leverages auxiliary information to enhance the statistical power of multiple hypothesis testing while controlling the false discovery rate (FDR) of high-dimensional data (such as 16S rRNA and WGS microbiome sequencing).

This repository contains all files related to the **chai** project.  
Please note that this repository does NOT include the `chai` R package itself.  
To install the `chai` package, please visit: **[chai](https://github.com/ziyiwang726/chai)**.

## Authors

**Ziyi Wang, Satabdi Saha, Christine B. Peterson, Yushu Shi**


## Repository Structure

- `data/`:
  Contains the real data used in the project.

- `scripts/`: 
  Contains all related R scripts and Python scripts.

- `reports/`:
  Contains the full workflow from data processing to visualization, including both code and generated figures.

- `figures/`:
  Contains all related figures.

## Setup Notes

### 1. Install the `chai` R package

The simulations and real-data scripts depend on the `chai` package, which is **not** included here. Install it from GitHub:

```r
devtools::install_github("ziyiwang726/chai")
```

The "main" Rmd files (e.g. `scripts/Realdata_Gastrectomy.Rmd`) call `library(chai)` and assume an installed copy. The "_YS" Rmd files (e.g. `scripts/Realdata_Gastrectomy_YS.Rmd`) also call `library(chai)` upstream; the only Rmd that still uses `devtools::load_all(...)` against a local source checkout is `scripts/JASA_Realdata_Melanoma_YS.Rmd`. To run that one, either:

- clone https://github.com/ziyiwang726/chai as a sibling of this project (so it lives at `../chai`), or
- export `CHAI_PKG_PATH=/path/to/chai` before launching R/RStudio.

### 2. Install the other R packages

```r
source("Gopalakrishnan/install_packages.R")
```

This script handles CRAN, Bioconductor, and GitHub-only dependencies (`CATMicrobiome`, `AdaPTGMM`, `adaptMT`, `FDRreg`, `OrderShapeEM`). Restart R after it finishes.

On macOS, install the official GNU Fortran toolchain from the R for macOS tools page before installing source-built packages such as `FDRreg`:
  https://mac.r-project.org/tools/

Current CRAN R builds expect the GNU Fortran runtime under `/opt/gfortran`. If it is missing, packages that link against Fortran libraries may fail during installation. The install script will detect this and print the exact `curl ... && sudo installer ...` command to run.

### 3. Working directory and data layout

All scripts expect to be run/knit **from the project root**. Data is laid out as:

```
data/Erawijantari/      # mtb.tsv, genera.counts.tsv, metadata.tsv
data/Zhu_autoencoder/   # ZhuF_2020_species_pathway.csv
data/Zhu_autoencoder/outputs/  # ZhuF_2020_species_AE5.csv, ZhuF_2020_species_PCA5.csv
data/chai_env/          # cached simulation .rds files (e.g. simA_all_results_with_zero.rds)
```

The melanoma 16S OTU table, taxonomy, metadata, and tree are loaded from the `CATMicrobiome` GitHub package via `system.file("extdata", ...)`, so no extra download is needed once that package is installed.

### 4. Known gaps in the YS scripts

The `_YS.Rmd` files are working drafts kept alongside the "main" Rmd files. A few sections still reference inputs that are not yet shipped in this repository — most notably:

- `scripts/JASA_Realdata_Melanoma_YS.Rmd` reads `./LLM_from_Yushu/d1Taxonomy_OTU_ABCD_class_and_below.csv`. The equivalent shipped file is `Gopalakrishnan/sideCov/combined/weighted/selectionTargetOTU/d1Taxonomy_OTU_ABCD_class_and_below_weighted.csv` (point the `read_csv(...)` call at it, or copy the file into `LLM_from_Yushu/`).

The "main" `Realdata_*.Rmd` files do not have this gap; they `load("./chai_env/<dataset>_env.RData")` to pick up pre-computed results, so you only need that `.RData` archive (not shipped — re-create it by knitting the YS variant once).

### 5. Source code layout

The R helper functions used by `chai` and by the analysis scripts live in `Gopalakrishnan/code/`. Earlier copies were duplicated under a top-level `code/` directory; that directory has now been removed and any incremental fixes (extra input validation, `+Inf` log-prob handling, NA/Inf checks for mixing weights) merged into the canonical `Gopalakrishnan/code/` versions.

### 6. LLM pipeline (Gopalakrishnan/LLMCode/) for external users

The literature-evidence pipeline under `Gopalakrishnan/LLMCode/` calls Anthropic / OpenAI / Gemini APIs. Copy `Gopalakrishnan/LLMCode/.env.example` to `Gopalakrishnan/LLMCode/.env` and fill in a key for whichever provider you plan to use; you do **not** need a Cornell key. See `Gopalakrishnan/LLMCode/README.md` for run commands.

## Simulation 1 - One dimensional auxiliary-information
For the generation and visualization of **Simulation 1**, please see: [Simulation 1](reports/JASA_Simulation_1.html).  

#### Simulation setting
Briefly, this simulation considers a multiple testing setting with 1,000 hypotheses, including 950 null hypotheses and two non-null groups with positive and negative signals. 

For each simulation run, $z$-statistics are generated from a three-component mixture: the null group follows $N(0,1)$, while the two alternative groups follow normal distributions centered at $2$ and $-2$ with variance 0.5. Two-sided $p$-values are then computed from the simulated $z$-statistics.

In Parallel, a side-information variable $\mathbf{X}$ is generated fro each hypothesis. The informativeness of $\mathbf{X}$ is controlled by a parameter $a$, where $a = 0$ means $\mathbf{X}$ is uninformative about signal status, and larger values of $a$ make the separation of $\mathbf{X}$ between null and non-null groups stronger. We considered $a \in {0, 0.5, 1, 1.5, 2}$, ranging from no information to highly informative side information. The simulation was repeated 100 times using different random seeds, and results were evaluated across a grid of target FDR levels from 0.01 to 0.10.

#### Results

Here is the line chart of **chai** and other benchmark methods in simulation 1 when informativeness parameter $a$ is fixed at 1 or 2:

![Simulation 1 - FDR and power](figures/SimA_fixa_1_2_zero.png)

## Simulation 2 - Two dimensional auxiliary-information
For the generation and visualization of **Simulation 2**, please see: [Simulation 2](reports/JASA_Simulation_2.html).  

#### Simulation setting
This simulation also considers a multiple testing setting with 1,000 hypotheses, including 950 null hypotheses and two non-null groups with positive (25) and negative signals (25). 

The null hypotheses were generated from standard normal distributions, while the two alternative groups were generated with positive and negative signals in both the test statistic and side infromation. Specifically, the side information included two variables, $x1$ and $x2$, which were combined into a two-dimensional covariate matrix $\mathbf{X}$. The simulation was repeated 100 times with different random seeds.


#### Results

This line chart shows the number of rejected hypotheses selected by **chai** and other benchmark methods in simulation 2 across different target FDR levels:

![Simulation 2 - FDR and power](figures/SimZ_fdr_power_updated_zero.png)

## Real data

### 1. Shotgun metagenomic sequencing data: gastrectomy vs. healthy individuals
For the full application and visualization of this real dataset, please see: [Gastrectomy](reports/JASA_Realdata_Gastrectomy.html).


### 2. 16S rRNA gene sequence data: responders vs. non-responders in a melanoma cohort
For the full application and visualization of this real dataset, please see: [Melanoma](reports/JASA_Realdata_Melanoma.html).


### 3. Shotgun metagenomic sequencing data: schizophrenia vs. healthy individuals
For the full application and visualization of this real dataset, please see: [Schizophrenia](reports/JASA_Realdata_Schizophrenia.html).
