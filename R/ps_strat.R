#' Propensity Score Stratification for Poisson Model
#'
#' Performs propensity score stratification and indirect standardization to
#' calculate observed and expected outcome counts for treated patients. This
#' method is used with TBSS Poisson probability models to adjust for confounding
#' by creating propensity score strata and calculating standardized expected counts.
#'
#' @param input_data A data frame containing patient data as returned by
#'   \code{plasmode()}. Must contain columns: \code{id}, \code{treatment_group},
#'   \code{outcome}, \code{ooi_label}, and covariate columns specified in the
#'   stratification formula.
#' @param strat_args A named list of arguments for stratification. Must include:
#'   \describe{
#'     \item{formula}{Formula specifying treatment ~ covariates for propensity
#'       score model (e.g., \code{treatment_group ~ age + sex})}
#'     \item{family}{GLM family, typically \code{binomial()} for binary treatment}
#'     \item{n_strata}{Integer number of propensity score strata to create (e.g., 5 or 10)}
#'   }
#'   Additional GLM arguments (e.g., \code{link}) can be included.
#'
#' @return A data frame with one row per outcome containing:
#'   \describe{
#'     \item{leaf}{Outcome codes}
#'     \item{expected}{Expected counts in treated group (standardized by referent group)}
#'     \item{observed}{Observed counts in treated group}
#'     \item{ooi_label}{Outcome of interest label}
#'   }
#'   This format is suitable for input to \code{TBSS::poissonTBSS()}.
#'
#' @details
#' **Stratification Process:**
#' \enumerate{
#'   \item Fits propensity score model using specified formula and family
#'   \item Predicts propensity scores for all patients
#'   \item Creates strata based on quantiles of fitted scores
#'   \item Counts outcomes by stratum in both treatment groups
#'   \item Calculates expected counts using indirect standardization:
#'     \deqn{E_k = \sum_s \frac{C_{0,s,k}}{N_{0,s}} \times N_{1,s}}
#'     where \eqn{C_{0,s,k}} is referent count for outcome k in stratum s,
#'     \eqn{N_{0,s}} is total referent in stratum s, and
#'     \eqn{N_{1,s}} is total treated in stratum s
#' }
#'
#' **Propensity Score Model:**
#' The propensity score is the probability of treatment given covariates,
#' estimated via logistic regression (or other GLM). Stratification on propensity
#' scores balances covariates between treatment groups within strata.
#'
#' **Indirect Standardization:**
#' Expected counts represent what would be observed in the treated group if it
#' had the same outcome distribution as the referent group, adjusted for
#' propensity score stratum composition.
#'
#' @export
ps_strat <- function(input_data, strat_args) {
  
  # Input validation
  
  # Check input_data is provided
  if (missing(input_data)) {
    stop("Argument 'input_data' is missing with no default.", call. = FALSE)
  }
  
  # Check strat_args is provided
  if (missing(strat_args)) {
    stop("Argument 'strat_args' is missing with no default.", call. = FALSE)
  }
  
  # Check input_data is a data frame
  if (!is.data.frame(input_data)) {
    stop("'input_data' must be a data frame or tibble.", call. = FALSE)
  }
  
  # Check input_data is not empty
  if (nrow(input_data) == 0) {
    stop("'input_data' is empty (contains no rows).", call. = FALSE)
  }
  
  # Check required columns exist
  required_cols <- c("id", "treatment_group", "outcome", "ooi_label")
  missing_cols <- setdiff(required_cols, names(input_data))
  
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "'input_data' is missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Check strat_args is a list
  if (!is.list(strat_args)) {
    stop("'strat_args' must be a list.", call. = FALSE)
  }
  
  # Check strat_args is not empty
  if (length(strat_args) == 0) {
    stop("'strat_args' cannot be empty.", call. = FALSE)
  }
  
  # Check required elements in strat_args
  required_args <- c("formula", "family", "n_strata")
  missing_args <- setdiff(required_args, names(strat_args))
  
  if (length(missing_args) > 0) {
    stop(
      sprintf(
        "'strat_args' is missing required element(s): %s",
        paste(missing_args, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Check formula is a formula object
  if (!inherits(strat_args$formula, "formula")) {
    stop("'strat_args$formula' must be a formula object.", call. = FALSE)
  }
  
  # Check n_strata is valid
  if (!is.numeric(strat_args$n_strata)) {
    stop("'strat_args$n_strata' must be numeric.", call. = FALSE)
  }
  
  if (length(strat_args$n_strata) != 1) {
    stop("'strat_args$n_strata' must be a single value.", call. = FALSE)
  }
  
  if (strat_args$n_strata < 2) {
    stop("'strat_args$n_strata' must be >= 2.", call. = FALSE)
  }
  
  n_strata <- as.integer(strat_args$n_strata)
  
  # Check both treatment groups exist
  unique_treatments <- unique(input_data$treatment_group)
  if (!all(c(0, 1) %in% unique_treatments)) {
    stop(
      sprintf(
        "'input_data' must contain both treatment groups (0 and 1).\nFound: %s",
        paste(unique_treatments, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Extract covariate names from formula
  formula_terms <- all.vars(strat_args$formula)
  covariate_cols <- setdiff(formula_terms, "treatment_group")
  
  # Check that covariates exist in data
  missing_covars <- setdiff(covariate_cols, names(input_data))
  if (length(missing_covars) > 0) {
    stop(
      sprintf(
        "Covariates specified in formula not found in 'input_data': %s\nAvailable columns: %s",
        paste(missing_covars, collapse = ", "),
        paste(names(input_data), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Prepare data for propensity score model
  
  # Extract unique patient records (one row per patient)
  cov <- input_data %>%
    dplyr::select(-outcome) %>%
    dplyr::distinct()
  
  # Convert to data frame (GLM requirement)
  cov_df <- as.data.frame(cov)
  
  # Check for missing values in covariates
  missing_in_covars <- sapply(cov_df[covariate_cols], function(x) sum(is.na(x)))
  if (any(missing_in_covars > 0)) {
    covars_with_na <- names(missing_in_covars[missing_in_covars > 0])
    warning(
      sprintf(
        "Covariates contain missing values:\n%s\nGLM may fail or drop observations.",
        paste(sprintf("  %s: %d missing", covars_with_na, missing_in_covars[covars_with_na]), 
              collapse = "\n")
      ),
      call. = FALSE
    )
  }
  
  # Fit propensity score model

  glm_args <- strat_args
  glm_args$n_strata <- NULL
  glm_args$data <- cov_df
  
  # Fit propensity score model
  ps_model <- tryCatch(
    {
      do.call(stats::glm, glm_args)
    },
    error = function(e) {
      stop(
        sprintf(
          "glm() failed with error:\n%s\n\nCheck your formula and family specification.",
          e$message
        ),
        call. = FALSE
      )
    }
  )
  
  # Check model convergence
  if (!ps_model$converged) {
    warning(
      "Propensity score model did not converge. Results may be unreliable.",
      call. = FALSE
    )
  }
  
  # Predict propensity scores
  
  fitted_scores <- tryCatch(
    {
      predict(ps_model, type = "response")
    },
    error = function(e) {
      stop(
        sprintf(
          "predict() failed with error:\n%s",
          e$message
        ),
        call. = FALSE
      )
    }
  )
  
  # Check for invalid propensity scores
  if (any(is.na(fitted_scores))) {
    warning(
      sprintf(
        "%d patients have NA propensity scores and will be excluded.",
        sum(is.na(fitted_scores))
      ),
      call. = FALSE
    )
  }
  
  if (any(fitted_scores < 0 | fitted_scores > 1, na.rm = TRUE)) {
    warning(
      "Some propensity scores are outside [0,1]. Check model specification.",
      call. = FALSE
    )
  }
  
  cov_df$fitted_score <- fitted_scores
  
  # Create propensity score strata
  
  # Calculate quantile breaks
  quantile_breaks <- quantile(
    fitted_scores,
    probs = seq(0, 1, 1 / n_strata),
    na.rm = TRUE
  )
  
  # Create strata
  cov_df$ps_stratum <- cut(
    cov_df$fitted_score,
    breaks = quantile_breaks,
    include.lowest = TRUE,
    labels = 1:n_strata
  )
  
  # Check stratum distribution
  stratum_counts <- table(cov_df$ps_stratum, cov_df$treatment_group)
  if (any(stratum_counts == 0)) {
    warning(
      "Some strata have zero patients in one or both treatment groups.\nConsider using fewer strata.",
      call. = FALSE
    )
  }
  
  # Calculate observed and expected counts
  
  # Get all unique outcomes
  unique_outcomes <- unique(input_data$outcome)
  
  # Create all outcome-stratum combinations
  all_combos <- tidyr::expand_grid(
    outcome = unique_outcomes,
    ps_stratum = factor(1:n_strata, levels = 1:n_strata)
  )
  
  # Join outcomes to patient data
  patient_outcomes <- cov_df %>%
    dplyr::select(id, treatment_group, ps_stratum) %>%
    dplyr::left_join(
      input_data %>% dplyr::select(id, outcome),
      by = "id",
      relationship = "many-to-many"
    )
  
  # Count referent group outcomes by stratum
  ref_counts <- patient_outcomes %>%
    dplyr::filter(treatment_group == 0) %>%
    dplyr::count(outcome, ps_stratum, name = "ref_count") %>%
    dplyr::right_join(all_combos, by = c("outcome", "ps_stratum")) %>%
    dplyr::mutate(ref_count = tidyr::replace_na(ref_count, 0))
  
  # Count treatment group outcomes by stratum
  treat_counts <- patient_outcomes %>%
    dplyr::filter(treatment_group == 1) %>%
    dplyr::count(outcome, ps_stratum, name = "treat_count") %>%
    dplyr::right_join(all_combos, by = c("outcome", "ps_stratum")) %>%
    dplyr::mutate(treat_count = tidyr::replace_na(treat_count, 0))
  
  # Calculate totals by stratum
  ref_totals <- ref_counts %>%
    dplyr::group_by(ps_stratum) %>%
    dplyr::mutate(total_ref = sum(ref_count)) %>%
    dplyr::ungroup()
  
  treat_totals <- treat_counts %>%
    dplyr::group_by(ps_stratum) %>%
    dplyr::mutate(total_treat = sum(treat_count)) %>%
    dplyr::ungroup()
  
  # Merge and calculate weighted expected counts
  merged <- ref_totals %>%
    dplyr::left_join(treat_totals, by = c("outcome", "ps_stratum")) %>%
    dplyr::mutate(
      treat_count = tidyr::replace_na(treat_count, 0),
      # Indirect standardization: (ref_rate * treat_total)
      weighted = dplyr::if_else(
        total_ref > 0,
        (ref_count / total_ref) * total_treat,
        0
      )
    )
  
  # Aggregate to outcome level
  result <- merged %>%
    dplyr::group_by(outcome) %>%
    dplyr::summarise(
      expected = round(sum(weighted, na.rm = TRUE)),
      observed = sum(treat_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(leaf = outcome)
  
  # Extract OOI label
  ooi_label <- cov_df$ooi_label[1]
  
  # Add OOI label to result
  result <- result %>%
    dplyr::mutate(ooi_label = ooi_label)
  
  return(result)
}

