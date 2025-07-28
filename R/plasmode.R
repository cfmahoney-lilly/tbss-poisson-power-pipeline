#' @title Plasmode simulation for randomly selected outcome.
#' @description Use plasmode simulation to enrich one outcome to the desired 
#'   level of relative risk.
#'   Specify exactly one of target_outcome (if you have pre-selected the outcome
#'   you wish to enrich) or target_inc (if you wish for the program to select a 
#'   random outcome close to the specified incidence rate).
#' @return A resampled dataset with a simulated outcome enriched to the desired 
#'   relative risk.
#' @param outcome_data Data with patient id, treatment status, covariates, and 
#'   outcomes, as returned by select_subsample().
#' @param target_outcome A user-specified outcome code.
#' @param target_inc A user-specified incidence rate.
#' @param rr Desired level of simulated relative risk.
plasmode <- function(outcome_data, target_outcome = NULL, target_inc = NULL,
                     rr = 2) {
  # select random ooi if necessary
  observed_incidence <- function(df, target_inc) {
    obs_inc <- df %>%
      count(outcome) %>%
      mutate(incidence = n / n_distinct(df$id)) %>%
      mutate(distance = abs(incidence - target_inc)) %>%
      slice_min(order_by = distance, n = 6) %>%
      dplyr::select(outcome, n, incidence) %>%
      arrange(desc(incidence))

    random_ooi <- obs_inc %>%
      dplyr::select(outcome) %>%
      slice_sample(n = 1) %>%
      pull()

    return(random_ooi)
  }

  if (is.null(target_outcome)) {
    target_outcome <- observed_incidence(outcome_data, target_inc)
  }

  # sample within treatment/comparator group
  boot_sample <- function(input_data, treatment) {
    n_id <- input_data %>%
      filter(treatment_group == treatment) %>%
      pull(id) %>%
      n_distinct()

    sample_cov <- input_data %>%
      filter(treatment_group == treatment) %>%
      group_by(id) %>%
      mutate(ooi = ifelse(any(outcome == target_outcome), 1, 0)) %>%
      ungroup() %>%
      dplyr::select(-outcome) %>%
      distinct() %>%
      slice_sample(n = n_id, replace = TRUE) %>%
      mutate(sim_id = row_number())

    return(sample_cov)
  }

  sample_df <- bind_rows(boot_sample(outcome_data, 1), boot_sample(outcome_data, 0)) %>%
    mutate(sim_id = row_number())

  # calculate probability of outcome given exposure and relative risk
  # exposed proportion in unmatched cohort
  exp_prob <- outcome_data %>%
    summarise(mean_treatment = mean(treatment_group)) %>%
    pull()

  # conditional probability of exposure given outcome
  binom_prob <- (rr * exp_prob) / ((rr * exp_prob) + exp_prob)

  # simulate exposure given that patient has ooi
  ooi_ids <- sample_df %>%
    filter(ooi == 1) %>%
    dplyr::select(sim_id, id) %>%
    mutate(exp_sim = rbinom(n = n(), size = 1, prob = binom_prob)) %>%
    mutate(perm_ids = sample(sim_id)) %>%
    dplyr::select(id, sim_id, exp_sim, perm_ids)

  ooi_sim_count <- ooi_ids %>%
    dplyr::select(sim_id, exp_sim) %>%
    distinct() %>%
    filter(exp_sim == 1) %>%
    summarise(count = n()) %>%
    pull(count)

  # remaining number to be assigned as exposed
  total_exp <- sample_df %>%
    filter(treatment_group == 1) %>%
    summarise(count = n()) %>%
    pull(count)

  non_ooi_exp <- total_exp - ooi_sim_count

  non_ooi_ids <- sample_df %>%
    filter(ooi == 0) %>%
    dplyr::select(id, sim_id) %>%
    distinct() %>%
    mutate(row_id = row_number()) %>%
    mutate(exp_sim = if_else(row_id %in% sample(row_id, non_ooi_exp), 1, 0)) %>%
    mutate(perm_ids = sample(sim_id)) %>%
    dplyr::select(id, sim_id, exp_sim, perm_ids)

  # simulated exposure vector
  sim_exp_df <- bind_rows(ooi_ids, non_ooi_ids)

  # join simulated exposure to incident outcomes
  sim_outcomes <- sim_exp_df %>%
    left_join(
      outcome_data %>%
        dplyr::select(id, treatment_group, outcome),
      by = "id",
      relationship = "many-to-many"
    ) %>%
    distinct()

  # join covariates and outcomes into single tibble
  # with ooi_label column

  cov <- sample_df %>%
    dplyr::select(-c(ooi, id)) %>%
    rename(id = sim_id)

  out <- sim_outcomes %>%
    dplyr::select(perm_ids, outcome) %>%
    rename(id = perm_ids)

  combined_df <- cov %>%
    left_join(out,
      by = "id"
    )

  final_df <- combined_df %>%
    mutate(ooi_label = target_outcome)

  return(final_df)
}
