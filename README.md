# Green public procurement and bidder participation

Replication materials for **“Environmental Criteria and Clauses in Public Procurement: Theory and Evidence on Their Effects on Bidder Participation”**, by Adrien Deschamps, François Maréchal and Pierre-Henri Morand.

The repository contains data and R scripts for the empirical analysis of environmental clauses, environmental award criteria and bidder participation in French public procurement in 2022–2023. Observations correspond to contract lots, or to contracts when they are not divided into lots. The data are derived from BOAMP procurement notices.

## Repository contents

| File | Description |
| --- | --- |
| `replication.R` | Main estimates, within-contract subsampling exercise, propensity-score reweighting analyses and robustness to controlling for award price. |
| `replication_additional.R` | Sector and buyer heterogeneity, buyer heterogeneity tests and robustness to relaxing the upper bid-count restriction. |
| `data_without_price.csv` | Baseline dataset, restricted to observations with ten offers or fewer. |
| `data_with_price.csv` | Dataset used for the specifications controlling for award price. |
| `data_without_price_extended.csv` | Extended dataset retaining additional observations with more than ten offers. |
| `README_ADDITIONAL.md` | Detailed documentation of the additional script and extended dataset. |
| `LICENSE` | License terms. |


## Software requirements

The original replication package specifies **R 4.5.1**. The scripts require the following packages:

```r
install.packages(c(
  "dplyr", "glmmTMB", "readr", "sandwich", "clubSandwich",
  "lmtest", "parameters", "nnet", "ggplot2", "future",
  "furrr", "performance"
))
```
