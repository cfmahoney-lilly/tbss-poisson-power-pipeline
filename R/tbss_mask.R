#' Tree-Based Scan Statistic (TBSS) Detection for Outcome of Interest
#'
#' Computes Tree-Based Scan Statistic using either Bernoulli or Poisson models
#' to detect statistical alerts in hierarchical outcome data. Returns a binary
#' indicator of whether the outcome of interest (OOI) is flagged as a statistical
#' alert at the specified significance level to preserve privacy in pharmacovigilance
#' applications.
#'
#' @param df A data frame containing matched or stratified patient data as returned
#'   by \code{ps_match()} (for Bernoulli) or \code{ps_strat()} (for Poisson).
#'   Must contain columns:
#'   \describe{
#'     \item{outcome}{Outcome codes}
#'     \item{treatment_group}{Treatment assignment (0 or 1)}
#'     \item{ooi_label}{Outcome of interest label}
#'     \item{exp_p}{Exposure proportion (required for Bernoulli model)}
#'   }
#' @param model Character string specifying the probability model. Must be either
#'   "bernoulli" or "poisson". Use "bernoulli" for binary detection (matched data)
#'   and "poisson" for count-based detection (stratified data). Default is "bernoulli".
#' @param tree_obj A data frame representing the hierarchical tree structure with
#'   columns \code{parent} and \code{child}, as returned by \code{tree_builder()}.
#'   Defines parent-child relationships between outcome codes.
#' @param sig Numeric value specifying the significance level (alpha) for detection.
#'   Must be between 0 and 1. Default is 0.05.
#'
#' @return Integer value (0 or 1):
#'   \describe{
#'     \item{1}{OOI detected as statistical alert (p-value < sig)}
#'     \item{0}{OOI not detected (p-value >= sig)}
#'   }
#'
#' @details
#' **TBSS Process:**
#' \enumerate{
#'   \item Formats data according to the specified model:
#'     \itemize{
#'       \item \strong{Bernoulli}: Direct matched group comparison
#'       \item \strong{Poisson}: Comparison on basis of underlying control rates
#'     }
#'   \item Computes tree-based scan statistic using the TBSS package
#'   \item Extracts results and filters by significance level
#'   \item Checks if the OOI appears in the significantly enriched nodes
#' }
#'
#' @export
tbss_mask <- function(df, model = "bernoulli", tree_obj, sig = 0.05) {
  
  # Input validation
  
  # Check df is provided
  if (missing(df)) {
    stop("Argument 'df' is missing with no default.", call. = FALSE)
  }
  
  # Check tree_obj is provided
  if (missing(tree_obj)) {
    stop("Argument 'tree_obj' is missing with no default.", call. = FALSE)
  }
  
  # Check df is a data frame
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame or tibble.", call. = FALSE)
  }
  
  # Check df is not empty
  if (nrow(df) == 0) {
    stop("'df' is empty (contains no rows).", call. = FALSE)
  }
  
  # Check model is valid
  if (!is.character(model) || length(model) != 1) {
    stop("'model' must be a single character string.", call. = FALSE)
  }
  
  model <- tolower(model)
  valid_models <- c("bernoulli", "poisson")
  
  if (!model %in% valid_models) {
    stop(
      sprintf(
        "'model' must be one of: %s\nProvided: '%s'",
        paste(valid_models, collapse = ", "),
        model
      ),
      call. = FALSE
    )
  }
  
  # Check tree_obj is a data frame
  if (!is.data.frame(tree_obj)) {
    stop("'tree_obj' must be a data frame or tibble.", call. = FALSE)
  }
  
  # Check tree_obj has required columns
  required_tree_cols <- c("parent", "child")
  missing_tree_cols <- setdiff(required_tree_cols, names(tree_obj))
  
  if (length(missing_tree_cols) > 0) {
    stop(
      sprintf(
        "'tree_obj' is missing required column(s): %s\nFound columns: %s",
        paste(missing_tree_cols, collapse = ", "),
        paste(names(tree_obj), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Check tree_obj is not empty
  if (nrow(tree_obj) == 0) {
    stop("'tree_obj' is empty (contains no rows).", call. = FALSE)
  }
  
  # Check sig is valid
  if (!is.numeric(sig)) {
    stop("'sig' must be numeric.", call. = FALSE)
  }
  
  if (length(sig) != 1) {
    stop("'sig' must be a single value.", call. = FALSE)
  }
  
  if (sig <= 0 || sig >= 1) {
    stop("'sig' must be between 0 and 1 (exclusive).", call. = FALSE)
  }
  
  # Check required columns in df
  required_cols <- c("outcome", "treatment_group", "ooi_label")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "'df' is missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Model-specific validation
  if (model == "bernoulli") {
    if (!"exp_p" %in% names(df)) {
      stop(
        "'df' must contain 'exp_p' column for Bernoulli model.",
        call. = FALSE
      )
    }
  }
  
  # Check TBSS package is available
  if (!requireNamespace("TBSS", quietly = TRUE)) {
    stop(
      "Package 'TBSS' is required but not installed.",
      call. = FALSE
    )
  }
  
  # Check both treatment groups exist
  unique_treatments <- unique(df$treatment_group)
  if (!all(c(0, 1) %in% unique_treatments)) {
    stop(
      sprintf(
        "'df' must contain both treatment groups (0 and 1).\nFound: %s",
        paste(unique_treatments, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Extract OOI (validate consistency)
  ooi_values <- unique(df$ooi_label)
  if (length(ooi_values) > 1) {
    warning(
      sprintf(
        "'ooi_label' column contains multiple values: %s\nUsing first value: '%s'",
        paste(ooi_values, collapse = ", "),
        ooi_values[1]
      ),
      call. = FALSE
    )
  }
  
  ooi <- ooi_values[1]
  
  # Check OOI exists in outcomes
  if (!ooi %in% df$outcome) {
    warning(
      sprintf(
        "Outcome of interest '%s' not found in outcome column.",
        ooi
      ),
      call. = FALSE
    )
  }
  
  # Run TBSS based on supplied model
  
  if (model == "bernoulli") {
    
    # Count outcomes by treatment group
    treat_outcomes <- df %>%
      dplyr::filter(treatment_group == 1) %>%
      dplyr::count(outcome, name = "case") %>%
      dplyr::rename(leaf = outcome)
    
    comp_outcomes <- df %>%
      dplyr::filter(treatment_group == 0) %>%
      dplyr::count(outcome, name = "control") %>%
      dplyr::rename(leaf = outcome)
    
    # Extract exposure proportion
    exp_p <- df$exp_p[1]
    
    # Validate exp_p
    if (is.na(exp_p)) {
      stop("'exp_p' column contains NA values.", call. = FALSE)
    }
    
    if (exp_p <= 0 || exp_p >= 1) {
      stop(
        sprintf(
          "'exp_p' must be between 0 and 1 (exclusive). Found: %.4f",
          exp_p
        ),
        call. = FALSE
      )
    }
    
    # Run Bernoulli TBSS
    mod <- tryCatch(
      {
        bernoulliTBSS(
          case = treat_outcomes,
          control = comp_outcomes,
          tree = tree_obj,
          p = exp_p
        )
      },
      error = function(e) {
        stop(
          sprintf(
            "bernoulliTBSS() failed with error:\n%s",
            e$message
          ),
          call. = FALSE
        )
      }
    )
    
  } else if (model == "poisson") {
    
    # Run Poisson TBSS
    mod <- tryCatch(
      {
        poissonTBSS(
          data = df_pois,
          tree = tree_obj
        )
      },
      error = function(e) {
        stop(
          sprintf(
            "poissonTBSS() failed with error:\n%s",
            e$message
          ),
          call. = FALSE
        )
      }
    )
  }
  
  # Extract results and check for OOI detection
  
  # Summarize and filter by significance
  mod_summary <- tryCatch(
    {
      summary(mod)
    },
    error = function(e) {
      stop(
        sprintf(
          "summary() on TBSS model failed with error:\n%s",
          e$message
        ),
        call. = FALSE
      )
    }
  )
  
  # Filter significant results
  sig_results <- mod_summary %>%
    dplyr::filter(pvalue < sig)
  
  # Check if OOI is reported in significant results
  if (ooi %in% sig_results$node) {
    return(1L)
  } else {
    return(0L)
  }
}
