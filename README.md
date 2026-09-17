# Target-trial-emulation-in-METEOR

This repository accompanies the manuscript:

**"Effect of biological DMARDs on physical function in patients with axial spondyloarthritis: results from a strategy trial emulation."**

The script provided here is the **analysis code used for the manuscript**. It produces the estimates, diagnostics, tables and figures reported in the article and its supplementary material.

## Contents

`Analysis in R.R` — a single annotated script, organised in the order in which results appear in the manuscript:

| Section | Output |
| --- | --- |
| Time-varying LTMLE, always versus never treated | Table 3; Supplementary Table S4 |
| Marginal structural model (`ipw`, `survey`) | Table 3; Supplementary Tables S5, S7, S8 |
| Parametric g-formula (`gfoRmula`) | Table 3; Supplementary Table S6 |
| Autoregressive time-lagged GEE (`geepack`) | Table 3 |
| LTMLE, all static and dynamic strategies, stratified on CRP | Figure 3; Supplementary Table S9 |
| Sensitivity analyses: censoring weights (IPCW) and alternative definitions of "treated" | Supplementary Tables S10 and S10.1 |

Table 1, Table 2 and Supplementary Table S3 were produced in Stata 15.1. Figure 1 and Supplementary Figure S1 were drawn with DAGitty and Excalidraw.

The script documents the data structure in detail: the variable-by-variable correspondence between the wide file used by LTMLE and the long file used by the MSM and the g-formula, the "next" convention by which a treatment variable at one visit refers to the interval that follows it, and which regressions each estimator actually fits.

## Data availability

The METEOR registry data cannot be shared, so **the script will not run as provided**. It is published so that the analytical choices behind the manuscript can be read, checked and reused, not so that the results can be re-executed.

Readers who wish to adapt the code to their own data will find the required column names and coding in the header of the script.

## Before running on your own data

* Replace `C:/mypath/` in the `setwd()` calls with your own directories. The script expects a `Datasets/` folder for input and `Tables/` and `Figures/` folders for output.
* The datasets read by the script are:
  * `meteor12ltmlemocondswitch.csv` — wide format, one row per patient, complete cases (n = 352)
  * `meteor12ltmlemocondswitchMI.csv` — wide format, all eligible patients (n = 761), used for the censoring-weighted analysis
  * `meteor12longnext.csv` — long format, three rows per patient, used by the MSM and the g-formula
  * `originalfulllong761.csv` — long format, supplies the alternative exposure definitions

## Software

Analyses were run in R 4.5.1 with:

| Package | Version | Used for |
| --- | --- | --- |
| `ltmle` | 1.3.0 | Longitudinal targeted maximum likelihood estimation |
| `ipw` | 1.2.1.1 | Numerator and denominator models for the treatment weights |
| `survey` | 4.5 | Marginal structural model with robust standard errors |
| `gfoRmula` | 1.1.1 | Parametric g-formula for time-varying treatments |
| `geepack` | 1.3.12 | Generalised estimating equations |

Also used: `data.table`, `dplyr`, `flextable`, `ggplot2`, `scales`, `writexl`.

Results are version-sensitive. In particular, the behaviour of `ltmle` around `gbounds` and `deterministic.g.function` should be checked against the version above before comparing output.

## Notes on two implementation choices

**Switching as an intervention.** In the three strategies that switch bDMARD, switching at 6 months is not a covariate but a third intervention node, placed between the two treatment decisions and given its own treatment model. Because only a patient who receives a bDMARD can switch, that model is constrained to a probability of zero whenever no bDMARD was received, using `deterministic.g.function`. Without the constraint the model assigns untreated patients a non-zero probability of switching.

**Weight truncation.** Cumulative probabilities are bounded, and the percentage of patients affected is reported as a positivity diagnostic. Throughout the manuscript this percentage refers to the strategy arm.

## Citation

If you use this repository, please cite the accompanying manuscript and the original methodological and software references, which are listed in the supplementary material.
