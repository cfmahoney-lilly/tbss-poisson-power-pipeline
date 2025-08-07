# Power calculation for TBSS with Poisson probability model
This pipeline uses the pipeline management package `targets` to streamline the evaluation of power in tree-based scan statistics used in pharmacovigilance. This particular pipeline is implemented using a Poisson probability model. A Bernoulli version is available in the repository [tbss-bernoulli-pipeline.](https://github.com/cfmahoney-lilly/tbss-bernoulli-power-pipeline)

## Background 
Tree-based scan statistics are a class of models that leverage hierarchical classification systems to identify increased risk of outcomes in an exposed group at multiple resolutions. They are used increasingly frequently in pharmacovigilance, but given their sensitivity to the tree-structured variable input and use of resampling-based multiple testing control, their statistical properties are complex to infer. In pharmacovigilance, it may be necessary to mask true patient data and potential statistical alerts, therefore this pipeline masks the output of the tree-based scan statistic and returns only a power value based on a binary indicator of whether a simulated outcome of interest alerted.

### Simulation Method 
The pipeline implements a plasmode simulation procedure based on that detailed in Wang, et al. (2018). Briefly, the outcome and covariate data are sampled with replacement, then a relative risk of known size is injected at the outcome of interest. Finally, treatment, outcome, and covariate data are permuted to preclude detection of unrelated signals present in the source data.

## Input data 
The specification of the precise study design, including inclusion/exclusion criteria, baseline period, and observation period are out of the scope of this pipeline. It is assumed that the user has created the cohorts according to their selected design and formatted the table as in *test_data.csv*, with a column for *id*, *treatment_group*, each relevant covariate, and a single *outcome* column with one row for each patient-outcome combination.  

| id | treatment_group | cov_1 | cov_2 | ... | outcome |
|----|-----------------|-------|-------|-----|---------|
| 1  | 1               | 0     | 1     | ... | A047    |
| 1  | 1               | 0     | 1     | ... | J459    |
| 1  | 1               | 0     | 1     | ... | E119    |
| 2  | 0               | 1     | 0     | ... | M545    |
| 2  | 0               | 1     | 0     | ... | F321    |


## Example implementation
For demonstration purposes, this repository contains sample simulated data and an ICD-10 sub-tree corresponding to all J-codes (respiratory system). The subset of codes is intended to reduce runtime in the example setting. The full curated ICD-10 tree as derived from the code list published by Sentinel is also available in the repository. 

The example pipeline is detailed in the _targets.R file, and each target and corresponding function or input contains comments to guide users in incorporating their own data and study design choices. Briefly, it selects treatment and comparator groups of the desired size, uses plasmode simulation to introduce known elevation in risk of an outcome, stratifies the data based on a user-defined propensity score function, and finally applies TBSS to determine if the simulated elevation in risk is detected. This process is repeated many times to determine the power, calculated as the proportion of simulations wherein the simulated excess risk is detected. The example uses only 10 simulations for demonstration; suitable numbers will vary by application, but 1,000 may be a good reference point.

The package `crew` allows the user to run computations in parallel. In the example, two workers are specified; the user can customise to their own computational resources.

Finally, after making all necessary modifications, the user calls `tar_make()` to run the pipeline.

## References 
Kulldorff, M., Dashevsky, I., Avery, T. R., Chan, A. K., Davis, R. L., Graham, D., ... & Brown, J. S. (2013). Drug safety data mining with a tree‐based scan statistic. Pharmacoepidemiology and drug safety, 22(5), 517-523.

Landau, W. M. (2021). The targets R package: a dynamic Make-like function-oriented pipeline toolkit for reproducibility and high-performance computing. Journal of Open Source Software, 6(57), 2959.

Russo, M., & Wang, S. V. (2024). An open‐source implementation of tree‐based scan statistics. Pharmacoepidemiology and Drug Safety, 33(3), e5765.

Wang, S. V., Maro, J. C., Baro, E., Izem, R., Dashevsky, I., Rogers, J. R., ... & Kulldorff, M. (2018). Data mining for adverse drug events with a propensity score-matched tree-based scan statistic. Epidemiology, 29(6), 895-903.
