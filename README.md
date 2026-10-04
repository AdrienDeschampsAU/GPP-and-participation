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

Keep the scripts and CSV files in the same repository folder. Both scripts can be run independently; `replication_additional.R` does not depend on objects created by `replication.R`.

## Software requirements

The original replication package specifies **R 4.5.1**. The scripts require the following packages:

```r
install.packages(c(
  "dplyr", "glmmTMB", "readr", "sandwich", "clubSandwich",
  "lmtest", "parameters", "nnet", "ggplot2", "future",
  "furrr", "performance"
))
```

The additional script alone requires only `dplyr`, `glmmTMB`, `readr`, `parameters` and `performance`. Package versions are not locked by this repository; retain the R and package versions used for each replication run.


frequency threshold is applied before splitting the baseline sample into sectors or buyer types.

The baseline and extended samples have the following checked sizes:

| Dataset | Rows in CSV | Rows after estimation filters | Rows above ten offers after filters |
| --- | ---: | ---: | ---: |
| `data_without_price.csv` | 52,510 | 49,214 | 0 |
| `data_without_price_extended.csv` | 53,529 | 50,214 | 1,000 |

The extended CSV uses the same variable names as the baseline dataset. It is a prepared estimation input, not a complete export of raw notices. The upstream contract-level exclusion described in the manuscript is not reconstructed by the additional script; reproducing that preparation step requires the complete data before filtering.

## Specifications and weighting

General specifications use zero-truncated negative-binomial models. Sector controls are either two-digit CPV fixed effects or three-digit CPV random intercepts. Quadratic specifications add `G_WEIGHT_SQUARED`.

Sampling weights are proportional to the inverse number of observations per contracting authority (`CAE_SIREN`) and normalized to sum to the sample size. Sector and buyer analyses recalculate these weights within each subsample.

- Sector heterogeneity uses the baseline sample, three-digit CPV **fixed effects**, and divisions 45 (construction), 90 (environmental services) and 33 (medical products).
- Buyer heterogeneity uses the baseline sample and three-digit CPV **random intercepts**. Buyer type is omitted as a control within a buyer-type subsample because it is constant there.
- Upper bid-count robustness uses the **entire extended sample**, including observations with ten offers or fewer. It does not estimate only on observations above ten offers.
- Propensity-score analyses combine treatment weights with sampling weights. The code distinguishes binary environmental instruments and categories of environmental-criterion weight.

Robust coefficient tables use `parameters::model_parameters(..., robust = TRUE)`, as implemented in the scripts. This call does not explicitly request contract-clustered standard errors. The buyer Q tests use the separate subgroup estimates and robust standard errors, assuming independence between subgroup coefficient estimates.

## Outputs

`replication.R` produces model output, subsampling summaries and a coefficient-distribution plot. When sourced in an R session, its fitted models and result objects remain available in that session; it does not provide a unified export directory.

`replication_additional.R` creates `results_additional/`:

| Output | Contents |
| --- | --- |
| `sample_audit.csv` | Observation counts at each cleaning stage. |
| `all_coefficients.csv` | Full robust coefficient tables with model and subgroup labels. |
| `green_coefficients.csv` | Green-procurement coefficients, standard errors, confidence intervals, p-values and exponentiated coefficients. |
| `model_diagnostics.csv` | Sample sizes, CPV groups, convergence status, Hessian checks and any errors. |
| `buyer_heterogeneity_Q_tests.csv` | Buyer heterogeneity tests for `G_CLAUSE` and `G_WEIGHT`. |
| `models.rds` | Fitted additional models. |
| `results_summary.txt` | Coefficient results and package-generated R-squared summaries. |
| `sessionInfo.txt` | R and package versions. |
