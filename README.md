# BRICS-EKC-Analysis

Replication materials for the study:

**“Calidad institucional, intensidad energética y sostenibilidad ambiental: Prueba de la hipótesis de la EKC en los BRICS bajo heterogeneidad política y estructural”**

This repository provides scripts to (i) download/assemble data (where licensing allows), (ii) build the BRICS panel (1990–2024), (iii) estimate EKC specifications under heterogeneity and cross-sectional dependence (e.g., CCE-MG / CS-ECM), and (iv) reproduce all tables and figures reported in the manuscript.

## Contents
- `scripts/` — data pipeline and estimation scripts
- `data/` — processed datasets (if shareable)
- `data_raw/` — raw downloads (not tracked)
- `tables/` — generated tables (outputs)
- `figures/` — generated figures (outputs)
- `docs/` — technical notes and supplementary material

## Requirements
- R (>= 4.3 recommended)
- Suggested packages: `tidyverse`, `plm`, `fixest`, `modelsummary`, `sandwich`, `lmtest`, `readr`, `WDI` (if World Bank API is used), and any package required by your scripts.

If you use `renv`, install dependencies with:
```r
install.packages("renv")
renv::restore()
