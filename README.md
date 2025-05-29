# 🐟 Fish Functional Rates in Mediterranean Sea

This project explores the use of **functions as rates** to describe key ecosystem processes in the Mediterranean Sea over the last 40 years, inspired by the idea that:
> **"A function is a rate."** [Jax 2005](https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/j.1600-0706.2005.13851.x)

We aim to define and compute three core ecological functions across fish species using data sourced primarily from FishBase and fishery landings. These functions are essential for quantifying ecosystem processes such as biomass production, trophic interactions, and nutrient cycling.
> The project is part of the [HORIZON B-USEFUL project](https://b-useful.eu).

---

## 📈 Core Functional Definitions

The project focuses on deriving species-level or group-level rates for the following ecological functions:

1. ### **Biomass Production (Growth Rate)**
   - Estimates the **rate of biomass accumulation** via species-specific growth parameters.
   - Typically based on von Bertalanffy growth function (VBGF) parameters (`L∞`, `K`, `t₀`) or derived growth performance indices.

2. ### **Predation/Consumption Rate**
   - Calculated as:  
     \[
     **Predation Rate** = (Q / B) × Biomass
     \]
   - Where:
     - `Q/B` = consumption-to-biomass ratio
     - `Biomass` = standing stock of the species or group

3. ### **Nutrient Cycling Rate**
   - Estimates the turnover of **Carbon (C), Nitrogen (N), and Phosphorus (P)** per unit of dry mass (DM), supplemented with **FishFlux** data (excretion rates).

---

## 📊 Data Sources

We intend to gather required parameters from the following sources:

- [FishBase](https://www.fishbase.se/)
  - Biological traits: growth rates, Q/B ratios, etc.
- [RFishBase](https://ropensci.github.io/rfishbase/)
  - R interface to FishBase (some data may not be directly accessible)
- [FishFlux](https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/1365-2435.13618)
  - Empirical estimates of nutrient excretion by fishes

> ⚠️ **Note:** Some key parameters may not be fully available through the RFishBase package. Manual data extraction, cross-referencing, or supplementary compilation may be necessary.

---

## 🛠️ Planned Workflow

1. **Define target species or communities**
2. **Extract growth and Q/B parameters**
3. **Estimate biomass production and consumption**
4. **Compile or model C:N:P excretion rates**
5. **Integrate and analyze functional rates across ecosystems**

---

## 📂 How to Navigate This Repository

This repository is organized to facilitate reproducibility and clarity in workflow. Here's how to find and use its contents:

- **📁 Data**  
  The data for this project will be released upon publication.  
  If you are part of the project team, you must manually add the datasets to the appropriate folder.  
  Please consult the `.gitignore` file to see which data files are expected but excluded from version control:  
  👉 [`.gitignore`](https://github.com/JayCrlt/Med_Fish_Functions/blob/main/.gitignore)

- **📊 Figures**  
  All generated figures, plots, and visual outputs will be stored here:  
  👉 [Outputs](https://github.com/JayCrlt/Med_Fish_Functions/tree/main/Outputs)

- **📜 Scripts**  
  Analysis and data processing scripts are stored and numbered sequentially in:  
  👉 [Scripts](https://github.com/JayCrlt/Med_Fish_Functions/tree/main/Scripts)

---

## 📦 Dependencies (R)

```r
─ Session info ────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.4.3 (2025-02-28)
 os       macOS Sequoia 15.5
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Madrid
 date     2025-05-29
 rstudio  2024.12.1+563 Kousa Dogwood (desktop)
 pandoc   NA
 quarto   1.5.57 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

─ Packages ────────────────────────────────────────────────────────────────────────────
 package           * version    date (UTC) lib source
 abind               1.4-8      2024-09-12 [1] CRAN (R 4.4.1)
 ade4                1.7-23     2025-02-14 [1] CRAN (R 4.4.1)
 ape               * 5.8-1      2024-12-16 [1] CRAN (R 4.4.1)
 arm                 1.14-4     2024-04-01 [1] CRAN (R 4.4.0)
 arrayhelpers        1.1-0      2020-02-04 [1] CRAN (R 4.4.0)
 backports           1.5.0      2024-05-23 [1] CRAN (R 4.4.1)
 blob                1.2.4      2023-03-17 [1] CRAN (R 4.4.0)
 boot                1.3-31     2024-08-28 [1] CRAN (R 4.4.3)
 broom               1.0.7      2024-09-26 [1] CRAN (R 4.4.1)
 cachem              1.1.0      2024-05-16 [1] CRAN (R 4.4.1)
 cellranger          1.1.0      2016-07-27 [1] CRAN (R 4.4.0)
 checkmate           2.3.2      2024-07-29 [1] CRAN (R 4.4.0)
 cli                 3.6.5      2025-04-23 [1] CRAN (R 4.4.1)
 clusterGeneration   1.3.8      2023-08-16 [1] CRAN (R 4.4.1)
 coda                0.19-4.1   2024-01-31 [1] CRAN (R 4.4.1)
 codetools           0.2-20     2024-03-31 [1] CRAN (R 4.4.3)
 colorspace          2.1-1      2024-07-26 [1] CRAN (R 4.4.1)
 combinat            0.0-8      2012-10-29 [1] CRAN (R 4.4.1)
 cowplot             1.1.3      2024-01-22 [1] CRAN (R 4.4.0)
 crayon              1.5.3      2024-06-20 [1] CRAN (R 4.4.1)
 curl                6.2.2      2025-03-24 [1] CRAN (R 4.4.1)
 DBI                 1.2.3      2024-06-02 [1] CRAN (R 4.4.1)
 dbplyr              2.5.0      2024-03-19 [1] CRAN (R 4.4.0)
 DEoptim             2.2-8      2022-11-11 [1] CRAN (R 4.4.1)
 Deriv               4.1.6      2024-09-13 [1] CRAN (R 4.4.1)
 deSolve             1.40       2023-11-27 [1] CRAN (R 4.4.1)
 digest              0.6.37     2024-08-19 [1] CRAN (R 4.4.1)
 distributional      0.5.0      2024-09-17 [1] CRAN (R 4.4.1)
 doBy                4.6.25     2025-01-30 [1] CRAN (R 4.4.1)
 doParallel          1.0.17     2022-02-07 [1] CRAN (R 4.4.0)
 dplyr             * 1.1.4      2023-11-17 [1] CRAN (R 4.4.0)
 duckdb              1.2.2      2025-04-29 [1] CRAN (R 4.4.1)
 duckdbfs            0.1.0      2025-04-04 [1] CRAN (R 4.4.1)
 expm                1.0-0      2024-08-19 [1] CRAN (R 4.4.1)
 farver              2.1.2      2024-05-13 [1] CRAN (R 4.4.1)
 fastmap             1.2.0      2024-05-15 [1] CRAN (R 4.4.1)
 fastmatch           1.1-6      2024-12-23 [1] CRAN (R 4.4.1)
 fishflux          * 0.0.1.7    2025-05-22 [1] Github (nschiett/fishflux@4d9d41e)
 fishtree          * 0.3.4      2021-01-31 [1] CRAN (R 4.4.0)
 fishualize          0.2.3      2022-03-08 [1] CRAN (R 4.4.0)
 forcats           * 1.0.0      2023-01-29 [1] CRAN (R 4.4.0)
 foreach             1.5.2      2022-02-02 [1] CRAN (R 4.4.0)
 fs                  1.6.6      2025-04-12 [1] CRAN (R 4.4.1)
 future              1.34.0     2024-07-29 [1] CRAN (R 4.4.0)
 future.apply        1.11.3     2024-10-27 [1] CRAN (R 4.4.1)
 geiger            * 2.0.11     2023-04-03 [1] CRAN (R 4.4.0)
 generics            0.1.4      2025-05-09 [1] CRAN (R 4.4.1)
 ggdist              3.3.3      2025-04-23 [1] CRAN (R 4.4.1)
 ggplot2           * 3.5.2      2025-04-09 [1] CRAN (R 4.4.1)
 ggridges          * 0.5.6      2024-01-23 [1] CRAN (R 4.4.0)
 globals             0.16.3     2024-03-08 [1] CRAN (R 4.4.1)
 glue                1.8.0      2024-09-30 [1] CRAN (R 4.4.1)
 gridExtra           2.3        2017-09-09 [1] CRAN (R 4.4.1)
 gtable              0.3.6      2024-10-25 [1] CRAN (R 4.4.1)
 hms                 1.1.3      2023-03-21 [1] CRAN (R 4.4.0)
 httr                1.4.7      2023-08-15 [1] CRAN (R 4.4.0)
 igraph              2.1.4      2025-01-23 [1] CRAN (R 4.4.1)
 inline              0.3.21     2025-01-09 [1] CRAN (R 4.4.1)
 iterators           1.0.14     2022-02-05 [1] CRAN (R 4.4.1)
 jsonlite            2.0.0      2025-03-27 [1] CRAN (R 4.4.1)
 labeling            0.4.3      2023-08-29 [1] CRAN (R 4.4.1)
 lattice             0.22-6     2024-03-20 [1] CRAN (R 4.4.3)
 lifecycle           1.0.4      2023-11-07 [1] CRAN (R 4.4.1)
 listenv             0.9.1      2024-01-29 [1] CRAN (R 4.4.1)
 lme4                1.1-36     2025-01-11 [1] CRAN (R 4.4.1)
 loo                 2.8.0      2024-07-03 [1] CRAN (R 4.4.0)
 lubridate         * 1.9.4      2024-12-08 [1] CRAN (R 4.4.1)
 magrittr            2.0.3      2022-03-30 [1] CRAN (R 4.4.1)
 maps              * 3.4.2.1    2024-11-10 [1] CRAN (R 4.4.1)
 MASS                7.3-64     2025-01-04 [1] CRAN (R 4.4.3)
 Matrix              1.7-2      2025-01-23 [1] CRAN (R 4.4.3)
 matrixStats         1.5.0      2025-01-07 [1] CRAN (R 4.4.1)
 memoise             2.0.1      2021-11-26 [1] CRAN (R 4.4.0)
 mgcv                1.9-1      2023-12-21 [1] CRAN (R 4.4.3)
 mi                  1.1        2022-06-06 [1] CRAN (R 4.4.0)
 microbenchmark      1.5.0      2024-09-04 [1] CRAN (R 4.4.1)
 minqa               1.2.8      2024-08-17 [1] CRAN (R 4.4.1)
 mnormt              2.1.1      2022-09-26 [1] CRAN (R 4.4.1)
 modelr              0.1.11     2023-03-22 [1] CRAN (R 4.4.0)
 mvtnorm             1.3-3      2025-01-10 [1] CRAN (R 4.4.1)
 nlme                3.1-167    2025-01-27 [1] CRAN (R 4.4.3)
 nloptr              2.1.1      2024-06-25 [1] CRAN (R 4.4.1)
 numDeriv            2016.8-1.1 2019-06-06 [1] CRAN (R 4.4.1)
 optimParallel       1.0-2      2021-02-11 [1] CRAN (R 4.4.1)
 parallelly          1.42.0     2025-01-30 [1] CRAN (R 4.4.1)
 patchwork         * 1.3.0      2024-09-16 [1] CRAN (R 4.4.1)
 phangorn            2.12.1     2024-09-17 [1] CRAN (R 4.4.1)
 phylobase           0.8.12     2024-01-30 [1] CRAN (R 4.4.0)
 phylolm             2.6.5      2024-09-30 [1] CRAN (R 4.4.1)
 phylopath           1.3.0      2024-06-11 [1] CRAN (R 4.4.0)
 phylosem          * 1.1.4      2024-04-02 [1] CRAN (R 4.4.0)
 phytools          * 2.4-4      2025-01-08 [1] CRAN (R 4.4.1)
 pillar              1.10.2     2025-04-05 [1] CRAN (R 4.4.1)
 pkgbuild            1.4.7      2025-03-24 [1] CRAN (R 4.4.3)
 pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.4.1)
 plyr                1.8.9      2023-10-02 [1] CRAN (R 4.4.1)
 png                 0.1-8      2022-11-29 [1] CRAN (R 4.4.1)
 posterior           1.6.1      2025-02-27 [1] CRAN (R 4.4.1)
 prettyunits         1.2.0      2023-09-24 [1] CRAN (R 4.4.1)
 progress            1.2.3      2023-12-06 [1] CRAN (R 4.4.0)
 purrr             * 1.0.4      2025-02-05 [1] CRAN (R 4.4.1)
 quadprog            1.5-8      2019-11-20 [1] CRAN (R 4.4.1)
 QuickJSR            1.7.0      2025-03-31 [1] CRAN (R 4.4.1)
 R6                  2.6.1      2025-02-15 [1] CRAN (R 4.4.1)
 ragg                1.3.3      2024-09-11 [1] CRAN (R 4.4.1)
 rbibutils           2.3        2024-10-04 [1] CRAN (R 4.4.1)
 RColorBrewer        1.1-3      2022-04-03 [1] CRAN (R 4.4.1)
 Rcpp                1.0.14     2025-01-12 [1] CRAN (R 4.4.1)
 RcppParallel        5.1.10     2025-01-24 [1] CRAN (R 4.4.1)
 Rdpack              2.6.2      2024-11-15 [1] CRAN (R 4.4.1)
 readr             * 2.1.5      2024-01-10 [1] CRAN (R 4.4.0)
 readxl            * 1.4.4      2025-02-27 [1] CRAN (R 4.4.1)
 reformulas          0.4.0      2024-11-03 [1] CRAN (R 4.4.1)
 reshape2            1.4.4      2020-04-09 [1] CRAN (R 4.4.0)
 rfishbase         * 5.0.1      2025-05-12 [1] Github (ropensci/rfishbase@eca5134)
 rlang               1.1.6      2025-04-11 [1] CRAN (R 4.4.1)
 rncl                0.8.7      2023-01-08 [1] CRAN (R 4.4.0)
 RNeXML              2.4.11     2023-02-01 [1] CRAN (R 4.4.0)
 Rphylopars        * 0.3.10     2024-01-18 [1] CRAN (R 4.4.0)
 rstan               2.32.7     2025-03-10 [1] CRAN (R 4.4.1)
 rstudioapi          0.17.1     2024-10-22 [1] CRAN (R 4.4.1)
 scales              1.4.0      2025-04-24 [1] CRAN (R 4.4.1)
 scatterplot3d       0.3-44     2023-05-05 [1] CRAN (R 4.4.1)
 sem                 3.1-16     2024-08-28 [1] CRAN (R 4.4.1)
 sessioninfo         1.2.3      2025-02-05 [1] CRAN (R 4.4.1)
 StanHeaders         2.32.10    2024-07-15 [1] CRAN (R 4.4.1)
 stringi             1.8.7      2025-03-27 [1] CRAN (R 4.4.1)
 stringr           * 1.5.1      2023-11-14 [1] CRAN (R 4.4.0)
 subplex             1.9        2024-07-05 [1] CRAN (R 4.4.1)
 svUnit              1.0.6      2021-04-19 [1] CRAN (R 4.4.1)
 systemfonts         1.2.1      2025-01-20 [1] CRAN (R 4.4.1)
 tensorA             0.36.2.1   2023-12-13 [1] CRAN (R 4.4.1)
 textshaping         1.0.0      2025-01-20 [1] CRAN (R 4.4.1)
 tibble            * 3.2.1      2023-03-20 [1] CRAN (R 4.4.0)
 tidybayes           3.0.7      2024-09-15 [1] CRAN (R 4.4.1)
 tidyr             * 1.3.1      2024-01-24 [1] CRAN (R 4.4.1)
 tidyselect          1.2.1      2024-03-11 [1] CRAN (R 4.4.0)
 tidyverse         * 2.0.0      2023-02-22 [1] CRAN (R 4.4.0)
 timechange          0.3.0      2024-01-18 [1] CRAN (R 4.4.1)
 TMB               * 1.9.17     2025-03-10 [1] CRAN (R 4.4.1)
 tzdb                0.5.0      2025-03-15 [1] CRAN (R 4.4.1)
 utf8                1.2.5      2025-05-01 [1] CRAN (R 4.4.1)
 uuid                1.2-1      2024-07-29 [1] CRAN (R 4.4.1)
 V8                  6.0.1      2025-02-02 [1] CRAN (R 4.4.1)
 vctrs               0.6.5      2023-12-01 [1] CRAN (R 4.4.0)
 withr               3.0.2      2024-10-28 [1] CRAN (R 4.4.1)
 XML                 3.99-0.18  2025-01-01 [1] CRAN (R 4.4.1)
 xml2                1.3.8      2025-03-14 [1] CRAN (R 4.4.1)

 [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
 * ── Packages attached to the search path.
 ```