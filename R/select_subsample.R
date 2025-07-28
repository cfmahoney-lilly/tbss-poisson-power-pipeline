#' @title Select random subsample of desired size for power calculation.
#' @description Select desired treated and comparator patient group sizes
#' @return A dataset with the selected number of treated and comparator patients
#' @param outcome_data Patient data as retured by get_data().
#' @param target_treated_n Number of treated patients to select.
#' @param target_comp_n Number of comparator patients to select.

select_subsample <- function(outcome_data, target_treated_n, target_comp_n) {
  sample_treat_ids <- outcome_data %>%
    filter(treatment_group == 1) %>%
    dplyr::select(id) %>%
    distinct() %>%
    slice_sample(n = target_treated_n) %>%
    pull()

  sample_treat <- outcome_data %>%
    filter(id %in% sample_treat_ids)

  sample_comp_ids <- outcome_data %>%
    filter(treatment_group == 0) %>%
    dplyr::select(id) %>%
    distinct() %>%
    slice_sample(n = target_comp_n) %>%
    pull()

  sample_comp <- outcome_data %>%
    filter(id %in% sample_comp_ids)

  sample_n <- bind_rows(sample_treat, sample_comp)

  return(sample_n)
}
