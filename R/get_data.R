#' @title Read in data in csv format
#' @description Data must be in csv format and contain an id variable called 
#'   "id", relevant covariate columns, a binary (treated == 1, untreated == 0)
#'   column indicating treatment status called "treatment_group", and a column 
#'   of outcome codes called "outcome"
#'   The table must have one row per patient-outcome, thus each row for a given
#'   patient id will have identical treatment_group and covariate values and a 
#'   unique outcome code.
#' @param file A text string of the name of the data file.

get_data <- function(file) {
  read_csv(file)
}
