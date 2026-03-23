# chai Project Repository
**chai** - Conditional Hypothesis testing using Auxiliary Information

**chai** is a covariate-informed statistical framework. It leverages auxiliary information to enhance the statistical power of multiple hypothesis testing while controlling the false discovery rate (FDR) of high-dimensional data (such as 16S rRNA and WGS microbiome sequencing).

This repository contains all files related to the **chai** project.  
Please note that this repository does NOT include the **chai** R package itself.  
To install the **chai** package, please visit: **[chai](https://github.com/ziyiwang726/chai)**.

## Authors

**Ziyi Wang, Satabdi Saha, Christine B. Peterson, Yushu Shi**


## Repository Structure

- `data/`:
  Contains the real data used in the project.

- `scripts/`: 
  Contains all related R scripts and Python scripts.

- `reports/`:
  Contains the full workflow from data processing to visualization, including both code and generated figures.

## Simulation 1 - One dimensional auxiliary-information
For the generation and visualization of **Simulation 1**, please see: [Simulation 1](reports/JASA_Simulation_1.html).  

Briefly, this simulation considers a multiple testing setting with 1,000 hypotheses, including 950 null hypotheses and two non-null groups with positive and negative signals. 

For each simulation run, $z$-statistics are generated from a three-component mixture: the null group follows $N(0,1)$, while the two alternative groups follow normal distributions centered at $2$ and $-2$ with variance 0.5. Two-sided $p$-values are then computed from the simulated $z$-statistics.

In Parallel, a side-information variable $\mathbf{X}$ is generated fro each hypothesis. The informativeness of $\mathbf{X}$ is controlled by a parameter $\alpha$, where $\alpha = 0$ means $\mathbf{X}$ is uninformative about signal status, and larger values of $\alpha$ make the separation between null and non-null groups stronger. We considered $\alpha \in {0, 0.5, 1, 1.5, 2}$, ranging from no information to highly informative side information. The simulation was repeated 100 times using different random seeds, and results were evaluated across a grid of target FDR levels from 0.01 to 0.10.


## Simulation 2 - Two dimensional auxiliary-information
For the generation and visualization of **Simulation 2**, please see: [Simulation 2](reports/JASA_Simulation_2.html).  
