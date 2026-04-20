# Gopalakrishnan

This folder contains scripts, helper functions, intermediate tables, and generated outputs used for the Gopalakrishnan microbiome and immunotherapy analysis workflow.

## Top-level scripts

- `install_prereqs_macos.sh`: macOS setup helper for system prerequisites.
- `install_packages.R`: installs the R packages needed by the analysis scripts.
- `real_data_16s_LLM_plot.R`: runs the main real-data 16S analysis using CHAI-related helper code and provider-specific `sideCov` inputs.
- `plot_class_family_order_LLM.R`: compares feature-selection behavior across taxonomic ranks and LLM/provider configurations.
- `plotDescriptivePlots.R`: builds descriptive plots from the `sideCov` result tables.

## Main folders

- `code/`: shared R helper functions sourced by the main analysis scripts, including CHAI utilities, conditional parameter handling, mixture-model code, and performance helpers.
- `LLMCode/`: pipeline code for literature retrieval, LLM-based extraction, matrix construction, and final table generation. This subfolder has its own README with more pipeline detail.
- `sideCov/`: provider-specific auxiliary result tables used as downstream inputs. The subfolders separate outputs from OpenAI, Claude, Gemini, combined runs, and older Gemini 3 runs.
- `cache/`: generated cache files created by the plotting and analysis scripts to avoid recomputation.
- `plots/`: generated figures and summary plot outputs.

## Notes

- Most scripts assume this folder is the working base directory and construct relative paths from the script location.
- `real_data_16s_LLM_plot.R` and `plot_class_family_order_LLM.R` both source helper functions from `code/`.
- Generated artifacts and local-only files such as `.env`, `.DS_Store`, `.Rhistory`, `cache/`, and `plots/` are excluded by this folder's `.gitignore`.
