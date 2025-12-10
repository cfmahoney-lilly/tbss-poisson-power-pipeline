#' @title Read in data in csv format
#' 
#' @description Reads and validates a CSV file containing patient outcome data.
#'   The function performs comprehensive validation to ensure data meets the 
#'   required structure for downstream analysis.
#'   
#' @param file A character string specifying the path to the CSV data file.
#' 
#' @return A tibble containing the validated data with the following structure:
#'   \describe{
#'     \item{id}{Unique patient identifier (can be any type)}
#'     \item{treatment_group}{Binary indicator (0 = untreated, 1 = treated)}
#'     \item{outcome}{Outcome codes (character or factor)}
#'     \item{...}{Additional covariate columns}
#'   }
#'   
#' @details 
#' The function validates that the input data meet the following requirements:
#' 
#' **File requirements:**
#' \itemize{
#'   \item File must be in CSV format (`.csv` extension)
#'   \item File must exist and be readable
#' }
#' 
#' **Required columns:**
#' \itemize{
#'   \item \code{id}: Patient identifier variable
#'   \item \code{treatment_group}: Binary treatment indicator (0 or 1)
#'   \item \code{outcome}: Outcome codes
#' }
#'
#' **Data validation:**
#' \itemize{
#'   \item \code{treatment_group} must contain only 0 (untreated) or 1 (treated)
#'   \item The file must contain at least one data row
#' }
#' 
#' @export
get_data <- function(file) {
  
  # Check file is provided
  if (missing(file)) {
    stop("Argument 'file' is missing with no default.", call. = FALSE)
  }
  
  # Check file is a character string
  if (!is.character(file)) {
    stop("'file' must be a character string (file path).", call. = FALSE)
  }
  
  # Check file is a single value
  if (length(file) != 1) {
    stop("'file' must be a single file path.", call. = FALSE)
  }
  
  # Check file is not empty string
  if (file == "") {
    stop("'file' cannot be an empty string.", call. = FALSE)
  }
  
  # Check file exists
  if (!file.exists(file)) {
    stop(sprintf("File not found: '%s'", file), call. = FALSE)
  }
  
  # Check file has .csv extension
  file_ext <- tolower(tools::file_ext(file))
  if (file_ext != "csv") {
    stop(
      sprintf(
        "File must be in CSV format. Found '.%s' extension.\nProvided file: '%s'",
        file_ext,
        file
      ),
      call. = FALSE
    )
  }
  
  # Check file is readable
  if (file.access(file, mode = 4) != 0) {
    stop(sprintf("File exists but is not readable: '%s'", file), call. = FALSE)
  }
  
  # Attempt to read CSV
  data <- tryCatch(
    {
      readr::read_csv(file, show_col_types = FALSE)
    },
    error = function(e) {
      stop(
        sprintf(
          "Error reading CSV file '%s': %s",
          file,
          e$message
        ),
        call. = FALSE
      )
    }
  )
  
  # Validate data structure
  
  # Check data is not empty
  if (nrow(data) == 0) {
    stop(sprintf("CSV file '%s' contains no data rows.", file), call. = FALSE)
  }
  
  # Check required columns exist
  required_cols <- c("id", "treatment_group", "outcome")
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "CSV file is missing required column(s): %s\nFound columns: %s",
        paste(missing_cols, collapse = ", "),
        paste(names(data), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Check treatment_group is binary (0 or 1)
  unique_treatments <- unique(data$treatment_group)
  unique_treatments <- unique_treatments[!is.na(unique_treatments)]
  
  if (!all(unique_treatments %in% c(0, 1))) {
    stop(
      sprintf(
        "'treatment_group' must be binary (0 = untreated, 1 = treated).\nFound values: %s",
        paste(sort(unique(data$treatment_group)), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Return validated data
  return(data)
}
