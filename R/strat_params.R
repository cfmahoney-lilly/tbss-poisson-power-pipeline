# Propensity score formula and stratification arguments
#   Additional arguments in the strat_args list will be passed to the glm()
#  This script contains only matching parameters and is not intended as a target
#  on its own.
#  Replace the propensity score formula with covariates relevant to your study
#  design; the sample covariates here correspond to categories of the Elixhauser
#  comorbidity score and are for representative purposes only

# propensity score formula
ps_form <- as.formula(treatment_group ~ gender + age + chf + carit + valv + pcd + 
                        pvd + hypunc + hypc + para + ond + cpd + diabunc + 
                        diabc + hypothy + rf + ld + pud + aids + lymph + 
                        metacanc + solidtum + rheumd + coag + obes + wloss + 
                        fed + blane + dane + alcohol + drug + psycho + depre)

# specify number of strata:
n_strata <- 10

#add glm parameters as required
strat_args <- list(
  formula = ps_form,
  family = binomial(),
  n_strata = n_strata
)
