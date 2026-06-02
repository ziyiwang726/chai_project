# chai Project Repository
 
**chai** — Conditional Hypothesis testing using Auxiliary Information
 
**chai** is a covariate-informed statistical framework that leverages auxiliary information to enhance the statistical power of multiple hypothesis testing while controlling the false discovery rate (FDR) in high-dimensional data (such as 16S rRNA and WGS microbiome sequencing).
 
> **Note:** This repository contains all analysis files related to the **chai** project, but does **not** include the `chai` R package itself.
> To install the package, please visit: [https://github.com/ziyiwang726/chai](https://github.com/ziyiwang726/chai)
 
---
 
## Authors
 
**Ziyi Wang, Satabdi Saha, Christine B. Peterson, Yushu Shi**
 
---
 
## Repository Structure
 
```
chai_project/
├── data/                   # Real datasets used in analyses
│   ├── Erawijantari/       # Gastrectomy metagenomic data
│   ├── Zhu_autoencoder/    # Schizophrenia WGS data + autoencoder outputs
│   └── chai_env/           # Cached simulation .rds files and .RData environments
├── scripts/                # R and Python analysis scripts (.Rmd, .R, .py)
├── reports/                # Full rendered workflows (HTML): code + figures
├── figures/                # Output figures
└── Gopalakrishnan/         # Melanoma 16S rRNA analysis pipeline
    ├── code/               # R helper functions used by analysis scripts
    ├── LLMCode/            # LLM-based literature evidence pipeline
    ├── cache/              # Cached intermediate results
    ├── plots/              # Intermediate plots and tables
    └── sideCov/            # Side covariate tables
```

---
 
## Setup
 
### 1. Install the `chai` R package
 
All simulation and real-data scripts depend on the `chai` package. Install it from GitHub:
 
```r
devtools::install_github("ziyiwang726/chai")
```

### 2. Install other R package dependencies
 
```r
source("scripts/install_packages.R")
```

This script handles CRAN, Bioconductor, and GitHub-only dependencies including `CATMicrobiome`, `AdaPTGMM`, `adaptMT`, `FDRreg`, and `OrderShapeEM`. **Restart R after it finishes.**

> **macOS users:** Install the official GNU Fortran toolchain from the [R for macOS tools page](https://mac.r-project.org/tools/) before installing source-built packages such as `FDRreg`. Current CRAN R builds expect the GNU Fortran runtime under `/opt/gfortran`. If it is missing, the install script will detect this and print the exact installation command to run.

### 3. Working directory and data layout
 
All scripts expect to be **run or knit from the project root** (`chai_project/`). Data is laid out as:
 
```
data/Erawijantari/              # mtb.tsv, genera.counts.tsv, metadata.tsv
data/Zhu_autoencoder/           # ZhuF_2020_species_pathway.csv
data/Zhu_autoencoder/outputs/   # ZhuF_2020_species_AE5.csv, ZhuF_2020_species_PCA5.csv
data/chai_env/                  # Cached simulation .rds files and dataset .RData environments
```
 
The melanoma 16S OTU table, taxonomy, metadata, and phylogenetic tree are loaded from the `CATMicrobiome` GitHub package via `system.file("extdata", ...)`, so no extra download is needed once that package is installed.
 
### 4. Pre-computed environments
 
The main `Realdata_*.Rmd` scripts load pre-computed results via `load("./chai_env/<dataset>_env.RData")` rather than re-running the full analysis. These `.RData` files are included in `data/chai_env/`. If you would like to fully re-generate them from scratch, you can either change the `eval` option in the chunks from `FALSE` to `TRUE` and then knit the file, or manually run the chunks.

### 5. LLM pipeline (`Gopalakrishnan/LLMCode/`)
 
The literature-evidence pipeline calls Anthropic, OpenAI, and Gemini APIs. To use it:
 
1. Copy `Gopalakrishnan/LLMCode/.env.example` to `Gopalakrishnan/LLMCode/.env`
2. Fill in an API key for whichever provider you plan to use (no Cornell key required)
3. See `Gopalakrishnan/LLMCode/README.md` for full run commands

---

## Simulation 1 - One dimensional auxiliary-information
For the generation and visualization of **Simulation 1**, please see: [Simulation 1](https://ziyiwang726.github.io/chai_project/reports/Simulation_1.html).  

#### Simulation setting
Briefly, this simulation considers a multiple testing setting with 1,000 hypotheses, including 950 null hypotheses and two non-null groups with positive and negative signals. 

For each simulation run, $z$-statistics are generated from a three-component mixture: the null group follows $N(0,1)$, while the two alternative groups follow normal distributions centered at $2$ and $-2$ with variance 0.5. Two-sided $p$-values are then computed from the simulated $z$-statistics.

In Parallel, a side-information variable $\mathbf{X}$ is generated fro each hypothesis. The informativeness of $\mathbf{X}$ is controlled by a parameter $a$, where $a = 0$ means $\mathbf{X}$ is uninformative about signal status, and larger values of $a$ make the separation of $\mathbf{X}$ between null and non-null groups stronger. We considered $a \in {0, 0.5, 1, 1.5, 2}$, ranging from no information to highly informative side information. The simulation was repeated 100 times using different random seeds, and results were evaluated across a grid of target FDR levels from 0.01 to 0.10.

#### Results

Here is the line chart of **chai** and other benchmark methods in simulation 1 when informativeness parameter $a$ is fixed at 1 or 2:

![Simulation 1 - FDR and power](figures/SimA_fixa_1_2_zero.png)

## Simulation 2 - Two dimensional auxiliary-information
For the generation and visualization of **Simulation 2**, please see: [Simulation 2](https://ziyiwang726.github.io/chai_project/reports/Simulation_2.html).  

#### Simulation setting
This simulation also considers a multiple testing setting with 1,000 hypotheses, including 950 null hypotheses and two non-null groups with positive (25) and negative signals (25). 

The null hypotheses were generated from standard normal distributions, while the two alternative groups were generated with positive and negative signals in both the test statistic and side infromation. Specifically, the side information included two variables, $x1$ and $x2$, which were combined into a two-dimensional covariate matrix $\mathbf{X}$. The simulation was repeated 100 times with different random seeds.


#### Results

This line chart shows the number of rejected hypotheses selected by **chai** and other benchmark methods in simulation 2 across different target FDR levels:

![Simulation 2 - FDR and power](figures/SimZ_fdr_power_updated_zero.png)

## Real data Analyses

### 1. Shotgun metagenomic sequencing data: gastrectomy vs. healthy individuals
Full workflow: [Gastrectomy](https://ziyiwang726.github.io/chai_project/reports/Realdata_Gastrectomy.html).

Canonical correlations between each genus and metabolite profiles are used as side information **X**.


### 2. 16S rRNA gene sequence data: responders vs. non-responders in a melanoma cohort
Full workflow: 
[Melanoma with PCoA](https://ziyiwang726.github.io/chai_project/reports/Realdata_Melanoma_PCoA.html) \\
[Melanoma (LLM)](https://ziyiwang726.github.io/chai_project/reports/Realdata_Melanoma_LLM.html)

Three analysis configurations were evaluated:
 
| Configuration | Test statistic | Side information **X** |
|---|---|---|
| (i) | Wilcoxon z-statistics | Phylogeny-derived PCoA coordinates | 
| (ii) | DESeq2 signed Wald statistics | Phylogeny-derived PCoA coordinates |
| (iii) | Wilcoxon z-statistics | LLM-derived covariates |


### 3. Shotgun metagenomic sequencing data: schizophrenia vs. healthy individuals
Full workflow: [Schizophrenia](https://ziyiwang726.github.io/chai_project/reports/Realdata_Schizophrenia.html).

Autoencoder-derived pathway covariates are used as side information **X**.
