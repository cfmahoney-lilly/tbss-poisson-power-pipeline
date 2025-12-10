#' Select Random Subsample for Power Calculation
#'
#' Randomly samples specified numbers of treated and comparator patients from
#' the dataset. All outcome rows for selected patients are retained to maintain
#' within-patient relationships.
#'
#' @param outcome_data A data frame or tibble containing patient outcome data
#'   as returned by \code{get_data()}. Must contain columns: \code{id},
#'   \code{treatment_group}, and \code{outcome}.
#' @param target_treated_n Integer specifying the number of treated patients
#'   (treatment_group == 1) to randomly sample. Must be a positive integer
#'   and cannot exceed the number of available treated patients.
#' @param target_comp_n Integer specifying the number of comparator patients
#'   (treatment_group == 0) to randomly sample. Must be a positive integer
#'   and cannot exceed the number of available comparator patients.
#'
#' @return A data frame containing all outcome rows for the randomly selected
#'   patients. The returned data maintains the same structure as the input,
#'   with one row per patient-outcome combination. 
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters
#'   \item Randomly samples \code{target_treated_n} unique patient IDs from
#'     the treated group (treatment_group == 1)
#'   \item Randomly samples \code{target_comp_n} unique patient IDs from
#'     the comparator group (treatment_group == 0)
#'   \item Returns all outcome rows for the selected patients from both groups
#' }
#'
#' The sampling is performed without replacement. Set a seed with \code{set.seed()}
#' before calling this function for reproducible results.
#'
#' @export
select_subsample <- function(outcome_data, target_treated_n, target_comp_n) {
  
  # Input validation
  
  # Check outcome_data is provided
  if (missing(outcome_data)) {
    stop("Argument 'outcome_data' is missing with no default.", call. = FALSE)
  }
  
  # Check target_treated_n is provided
  if (missing(target_treated_n)) {
    stop("Argument 'target_treated_n' is missing with no default.", call. = FALSE)
  }
  
  # Check target_comp_n is provided
  if (missing(target_comp_n)) {
    stop("Argument 'target_comp_n' is missing with no default.", call. = FALSE)
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
  required_cols <- c("id", "treatment_group")
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
  
  # Validate target_treated_n
  if (!is.numeric(target_treated_n)) {
    stop("'target_treated_n' must be numeric.", call. = FALSE)
  }
  
  if (length(target_treated_n) != 1) {
    stop("'target_treated_n' must be a single value.", call. = FALSE)
  }
  
  if (target_treated_n < 1) {
    stop("'target_treated_n' must be a positive integer (>= 1).", call. = FALSE)
  }
  
  target_treated_n <- as.integer(target_treated_n)
  
  # Validate target_comp_n
  if (!is.numeric(target_comp_n)) {
    stop("'target_comp_n' must be numeric.", call. = FALSE)
  }
  
  if (length(target_comp_n) != 1) {
    stop("'target_comp_n' must be a single value.", call. = FALSE)
  }
  
  if (target_comp_n < 1) {
    stop("'target_comp_n' must be a positive integer (>= 1).", call. = FALSE)
  }
  
  target_comp_n <- as.integer(target_comp_n)
  
  # Check available sample sizes 
  
  # Get unique patient counts by treatment group
  patient_counts <- outcome_data %>%
    dplyr::distinct(id, treatment_group) %>%
    dplyr::count(treatment_group, name = "n_patients")
  
  n_treated <- patient_counts %>%
    dplyr::filter(treatment_group == 1) %>%
    dplyr::pull(n_patients)
  
  n_comp <- patient_counts %>%
    dplyr::filter(treatment_group == 0) %>%
    dplyr::pull(n_patients)
  
  # Handle case where treatment groups have zero patients
  if (length(n_treated) == 0) {
    n_treated <- 0
  }
  
  if (length(n_comp) == 0) {
    n_comp <- 0
  }
  
  # Verify that sufficient treated patients are available
  if (n_treated < target_treated_n) {
    stop(
      sprintf(
        "Insufficient treated patients available.\nRequested: %d, Available: %d",
        target_treated_n,
        n_treated
      ),
      call. = FALSE
    )
  }
  
  # Verify that sufficient comparator patients are available
  if (n_comp < target_comp_n) {
    stop(
      sprintf(
        "Insufficient comparator patients available.\nRequested: %d, Available: %d",
        target_comp_n,
        n_comp
      ),
      call. = FALSE
    )
  }
  
  # Random smapling
  
  # Unique patient IDs by treatment group
  unique_patients <- outcome_data %>%
    dplyr::distinct(id, treatment_group)
  
  # Sample treated patient IDs
  sample_treat_ids <- unique_patients %>%
    dplyr::filter(treatment_group == 1) %>%
    dplyr::slice_sample(n = target_treated_n) %>%
    dplyr::pull(id)
  
  # Sample comparator patient IDs
  sample_comp_ids <- unique_patients %>%
    dplyr::filter(treatment_group == 0) %>%
    dplyr::slice_sample(n = target_comp_n) %>%
    dplyr::pull(id)
  
  # Combine sampled IDs
  selected_ids <- c(sample_treat_ids, sample_comp_ids)
  
  # Filter original data using combined IDs 
  sample_data <- outcome_data %>%
    dplyr::filter(id %in% selected_ids)
  
  return(sample_data)
}

