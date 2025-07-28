#' @title Propensity score stratification for Poisson probability model.
#' @description Propensity score sratification and indirect stratification for 
#'   use in TBSS with Poisson probability model. Arguments are specified in the 
#'   file strat_params.R
#' @return Observed and expected outcome counts for treated patients
#' @param input_data Data with patient id, treatment status, covariates, and 
#'   outcomes, as returned by plasmode().
#' @param match_args List of arguments for matching algorithm as specified in 
#'   strat_params.R.
ps_strat <- function(input_data, strat_args) {
  cov <- input_data %>%
    dplyr::select(-outcome) %>%
    distinct()

  cov_df <- as.data.frame(cov)

  glm_args <- strat_args
  glm_args$n_strata <- NULL

  args_df <- c(glm_args, list(data = cov_df))

  browser()
  # fit glm with specified formula
  ps_model <- do.call(glm, args_df)
  cov_df$fitted_score <- predict(ps_model, type = "response")

  cov_df <- cov_df %>% # stratify referent group
    mutate(ps_stratum = cut(fitted_score,
      breaks = quantile(fitted_score,
        probs = seq(0, 1, 1 / strat_args$n_strata)
      ),
      include.lowest = TRUE,
      labels = 1:strat_args$n_strata
    ))


  all_combos <- expand.grid( # accommodate 0 values
    outcome = unique(input_data$outcome),
    ps_stratum = factor(1:strat_args$n_strata,
      levels = 1:strat_args$n_strata
    )
  )

  # counts per outcome and stratum in referent group
  ref_counts <- cov_df %>%
    filter(treatment_group == 0) %>%
    left_join(input_data %>% select(id, outcome), by = "id") %>%
    count(outcome, ps_stratum, name = "ref_count") %>%
    right_join(all_combos, by = c("outcome", "ps_stratum")) %>%
    mutate(ref_count = replace_na(ref_count, 0))

  # counts in treatment group
  treat_counts <- cov_df %>%
    filter(treatment_group == 1) %>%
    left_join(input_data %>% select(id, outcome), by = "id") %>%
    count(outcome, ps_stratum, name = "treat_count") %>%
    right_join(all_combos, by = c("outcome", "ps_stratum")) %>%
    mutate(treat_count = replace_na(treat_count, 0))

  ref_totals <- ref_counts %>%
    group_by(ps_stratum) %>%
    mutate(total_ref = sum(ref_count)) %>%
    ungroup()

  treat_totals <- treat_counts %>%
    group_by(ps_stratum) %>%
    mutate(total_treat = sum(treat_count)) %>%
    ungroup()

  merged <- ref_totals %>%
    left_join(treat_totals, by = c("outcome", "ps_stratum")) %>%
    mutate(treat_count = replace_na(treat_count, 0)) %>%
    mutate(weighted = (ref_count / total_ref) * total_treat)

  result <- merged %>%
    group_by(outcome) %>%
    summarise(
      expected = round(sum(weighted, na.rm = TRUE)),
      observed = sum(treat_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename("leaf" = "outcome") %>%
    mutate(ooi_label = cov_df %>%
      slice(1) %>%
      pull(ooi_label))

  return(result)
}
