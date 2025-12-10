#' Plasmode Simulation for Outcome Enrichment
#'
#' Performs plasmode simulation to enrich a specified outcome (or randomly
#' selected outcome) to a desired relative risk level through bootstrap
#' resampling while preserving the relationship between covariates and outcomes.
#'
#' @param outcome_data A data frame or tibble containing patient outcome data
#'   as returned by \code{select_subsample()} or \code{get_data()}. Must contain
#'   columns: \code{id}, \code{treatment_group}, and \code{outcome}.
#' @param target_outcome Character string specifying the outcome code to enrich.
#'   If NULL, the function will randomly select an outcome near \code{target_inc}.
#'   Exactly one of \code{target_outcome} or \code{target_inc} must be specified.
#'   Default is NULL.
#' @param target_inc Numeric value specifying the desired incidence rate
#'   (between 0 and 1) for outcome selection. The function will randomly select
#'   one of the 6 outcomes closest to this incidence rate. If NULL, 
#'   \code{target_outcome} must be specified. Default is NULL.
#' @param rr Numeric value specifying the desired relative risk for the enriched
#'   outcome. Must be positive. Values > 1 indicate increased risk in the treated
#'   group; values < 1 indicate decreased risk. Default is 2.
#'
#' @return A data frame with the following structure:
#'   \describe{
#'     \item{id}{Simulated patient identifiers (sequential integers)}
#'     \item{treatment_group}{Simulated treatment assignment (0 or 1)}
#'     \item{...}{All covariate columns from input data}
#'     \item{outcome}{Outcome codes from resampled patients}
#'     \item{ooi_label}{Character column containing the enriched outcome code}
#'   }
#'   The returned data maintains one row per patient-outcome, with the outcome
#'   of interest enriched to achieve the specified relative risk.
#'
#' @details
#' As a sketch, the plasmode simulation method:
#' \enumerate{
#'   \item Selects or identifies the outcome of interest (OOI)
#'   \item Bootstrap resamples patients within treatment groups (with replacement)
#'   \item Calculates probability of treatment given outcome based on desired RR
#'   \item Simulates treatment assignment to achieve target relative risk
#'   \item Permutes patient IDs to preserve positivity
#'   \item Returns enriched dataset with preserved covariate-outcome relationships
#' }
#'
#' **Outcome Selection:**
#' \itemize{
#'   \item If \code{target_outcome} is specified, that outcome is enriched
#'   \item If \code{target_inc} is specified, the function finds the 6 outcomes
#'     with incidence closest to the target, then randomly selects one
#'   \item Incidence is calculated as: (patients with outcome) / (total patients)
#' }
#'
#' **Relative Risk Calculation:**
#' Uses the formula: P(Exposure|Outcome) = (RR × P(Exposure)) / (RR × P(Exposure) + (1 - P(Exposure)))
#'
#' @export
plasmode <- function(outcome_data, target_outcome = NULL, target_inc = NULL, rr = 2) {
  
  # Input validation
  
  # Check outcome_data is provided
  if (missing(outcome_data)) {
    stop("Argument 'outcome_data' is missing with no default.", call. = FALSE)
  }
  
  # Check outcome_data is a data frame
  if (!is.data.frame(outcome_data)) {
    stop("'outcome_data' must be a data frame or tibble.", call. = FALSE)
  }
  
  # Check outcome_data is not empty
  if (nrow(outcome_data) == 0) {
    stop("'outcome_data' is empty (contains no rows).", call. = FALSE)
  }
  
  # Check required columns exist
  required_cols <- c("id", "treatment_group", "outcome")
  missing_cols <- setdiff(required_cols, names(outcome_data))
  
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "'outcome_data' is missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Validate exactly one of target_outcome or target_inc is specified
  if (is.null(target_outcome) && is.null(target_inc)) {
    stop(
      "Exactly one of 'target_outcome' or 'target_inc' must be specified.",
      call. = FALSE
    )
  }
  
  if (!is.null(target_outcome) && !is.null(target_inc)) {
    stop(
      "Cannot specify both 'target_outcome' and 'target_inc'. Select one.",
      call. = FALSE
    )
  }
  
  # Validate target_outcome if provided
  if (!is.null(target_outcome)) {
    
    # Check target_outcome exists in data
    if (!target_outcome %in% outcome_data$outcome) {
      stop(
        sprintf(
          "'target_outcome' value '%s' not found in outcome_data.\nAvailable outcomes: %s",
          target_outcome,
          paste(head(unique(outcome_data$outcome), 10), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  
  # Validate target_inc if provided
  if (!is.null(target_inc)) {
    if (!is.numeric(target_inc)) {
      stop("'target_inc' must be numeric.", call. = FALSE)
    }
    
    if (length(target_inc) != 1) {
      stop("'target_inc' must be a single value.", call. = FALSE)
    }
    
    if (target_inc <= 0 || target_inc >= 1) {
      stop("'target_inc' must be between 0 and 1 (exclusive).", call. = FALSE)
    }
  }
  
  # Validate rr
  if (!is.numeric(rr)) {
    stop("'rr' must be numeric.", call. = FALSE)
  }
  
  if (length(rr) != 1) {
    stop("'rr' must be a single value.", call. = FALSE)
  }
  
  if (rr <= 0) {
    stop("'rr' must be positive (> 0).", call. = FALSE)
  }
  
  # Check both treatment groups exist
  unique_treatments <- unique(outcome_data$treatment_group)
  if (!all(c(0, 1) %in% unique_treatments)) {
    stop(
      sprintf(
        "'outcome_data' must contain both treatment groups (0 and 1).\nFound: %s",
        paste(unique_treatments, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Select outcome of interest
  
  if (is.null(target_outcome)) {
    # Calculate incidence for all outcomes
    n_patients <- dplyr::n_distinct(outcome_data$id)
    
    obs_inc <- outcome_data %>%
      dplyr::count(outcome, name = "n") %>%
      dplyr::mutate(
        incidence = n / n_patients,
        distance = abs(incidence - target_inc)
      ) %>%
      dplyr::slice_min(order_by = distance, n = 6) %>%
      dplyr::arrange(dplyr::desc(incidence))
    
    # Randomly select one of the 6 closest values
    target_outcome <- obs_inc %>%
      dplyr::slice_sample(n = 1) %>%
      dplyr::pull(outcome)
    
    message(sprintf(
      "Selected outcome: '%s' (incidence: %.4f, target: %.4f)",
      target_outcome,
      obs_inc %>% dplyr::filter(outcome == target_outcome) %>% dplyr::pull(incidence),
      target_inc
    ))
  }
  
  # Calculate patient-level data
  patient_data <- outcome_data %>%
    dplyr::group_by(id, treatment_group) %>%
    dplyr::summarise(
      ooi = as.integer(any(outcome == target_outcome)),
      .groups = "drop"
    )
  
  # Join all covariates 
  covariate_cols <- setdiff(names(outcome_data), c("outcome"))
  
  patient_data <- outcome_data %>%
    dplyr::select(dplyr::all_of(covariate_cols)) %>%
    dplyr::distinct() %>%
    dplyr::left_join(patient_data, by = c("id", "treatment_group"))
  
  # Bootstrap sample within treatment groups
  
  # Sample treated patients
  n_treated <- sum(patient_data$treatment_group == 1)
  sample_treated <- patient_data %>%
    dplyr::filter(treatment_group == 1) %>%
    dplyr::slice_sample(n = n_treated, replace = TRUE)
  
  # Sample comparator patients
  n_comp <- sum(patient_data$treatment_group == 0)
  sample_comp <- patient_data %>%
    dplyr::filter(treatment_group == 0) %>%
    dplyr::slice_sample(n = n_comp, replace = TRUE)
  
  # Combine and assign simulation IDs
  sample_df <- dplyr::bind_rows(sample_treated, sample_comp) %>%
    dplyr::mutate(sim_id = dplyr::row_number())
  
  # Calculate probability of exposure given outcome
  
  # Proportion exposed in original cohort
  exp_prob <- mean(outcome_data$treatment_group)
  
  # Conditional probability: P(Exposure|Outcome)
  binom_prob <- (rr * exp_prob) / ((rr * exp_prob) + (1 - exp_prob))
  
  # Simulate treatment assignment ----
  
  # For patients with OOI, simulate exposure based on RR
  ooi_patients <- sample_df %>%
    dplyr::filter(ooi == 1) %>%
    dplyr::select(id, sim_id)
  
  n_ooi <- nrow(ooi_patients)
  
  ooi_ids <- ooi_patients %>%
    dplyr::mutate(
      exp_sim = rbinom(n = n_ooi, size = 1, prob = binom_prob),
      perm_ids = sample(sim_id)
    )
  
  # Count how many OOI patients assigned to treatment
  ooi_sim_count <- sum(ooi_ids$exp_sim)
  
  # Remaining exposed slots to fill from non-OOI patients
  total_exp <- n_treated
  non_ooi_exp <- total_exp - ooi_sim_count
  
  # Account for edge case where all treatment designations filled by OOI patients
  if (non_ooi_exp < 0) {
    
    # Randomly deselect some OOI patients for treatment
    n_deselect <- abs(non_ooi_exp)
    ooi_ids <- ooi_ids %>%
      dplyr::mutate(
        row_id = dplyr::row_number(),
        exp_sim = dplyr::if_else(
          exp_sim == 1 & row_id %in% sample(which(exp_sim == 1), n_deselect),
          0L,
          exp_sim
        )
      ) %>%
      dplyr::select(-row_id)
    
    non_ooi_exp <- 0
  }
  
  # Assign exposure to non-OOI patients
  non_ooi_patients <- sample_df %>%
    dplyr::filter(ooi == 0) %>%
    dplyr::select(id, sim_id)
  
  if (nrow(non_ooi_patients) > 0 && non_ooi_exp > 0) {
    # Randomly select non-OOI patients for treatment
    selected_indices <- sample(seq_len(nrow(non_ooi_patients)), non_ooi_exp)
    
    non_ooi_ids <- non_ooi_patients %>%
      dplyr::mutate(
        row_id = dplyr::row_number(),
        exp_sim = as.integer(row_id %in% selected_indices),
        perm_ids = sample(sim_id)
      ) %>%
      dplyr::select(id, sim_id, exp_sim, perm_ids)
  } else if (nrow(non_ooi_patients) > 0) {
    # No additional exposed needed
    non_ooi_ids <- non_ooi_patients %>%
      dplyr::mutate(
        exp_sim = 0L,
        perm_ids = sample(sim_id)
      ) %>%
      dplyr::select(id, sim_id, exp_sim, perm_ids)
  } else {
    # Account for edge case of 0 non-OOI patients 
    non_ooi_ids <- ooi_ids %>%
      dplyr::slice(0) %>%
      dplyr::select(id, sim_id, exp_sim, perm_ids)
  }
  
  # Combine simulated exposure assignments
  sim_exp_df <- dplyr::bind_rows(ooi_ids, non_ooi_ids) %>%
    dplyr::select(id, sim_id, exp_sim, perm_ids)
  
  # Join simulated exposure to outcomes 
  
  # Create outcome lookup
  outcome_lookup <- outcome_data %>%
    dplyr::select(id, outcome) %>%
    dplyr::distinct()
  
  sim_outcomes <- sim_exp_df %>%
    dplyr::left_join(outcome_lookup, by = "id", relationship = "many-to-many") %>%
    dplyr::distinct()
  
  # Construct final output dataset
  
  # Covariates with simulated treatment assignment
  cov <- sample_df %>%
    dplyr::select(-ooi, -id) %>%
    dplyr::rename(id = sim_id)
  
  # Outcomes with permuted IDs
  out <- sim_outcomes %>%
    dplyr::select(perm_ids, outcome) %>%
    dplyr::rename(id = perm_ids)
  
  # Combine and add OOI label
  final_df <- cov %>%
    dplyr::left_join(out, by = "id") %>%
    dplyr::mutate(ooi_label = target_outcome)
  
  return(final_df)
}

