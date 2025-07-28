#' @title Format code list as hierarchical tree.
#' @description Accepts a vector of codes and generates a hierarchical tree with
#'  columns "child" and "parent"
#'  It is not necessary to remove nodes above leaf level. The function will 
#'    automatically de-duplicate parent-child relationships
#' @return A data frame of a hierarchical tree suitable for use with TBSS
#' @param code_list A vector of codes from a hierarchical coding system 
#' @param min_level The highest level for which to produce a parent node, for
#'  example, ICD-10 level 3 (A00) 
tree_builder <- function(code_list, min_level) {
  require(dplyr)
  require(stringr)
  
  tree <- tibble(child = character(), parent = character())
  
  for (code in code_list) {
    while (nchar(code) > min_level) {
      parent_code <- str_sub(code, 1, -2)
      tree <- tree %>% add_row(
        child = code,
        parent = parent_code
      )
      
      code <- parent_code
    }
  }
  
  tree_unique <- tree %>% distinct()
  
  return(tree_unique)
}
