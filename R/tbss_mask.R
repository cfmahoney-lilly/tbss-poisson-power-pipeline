#' @title TBSS binary detection flag for specified outcome.
#' @description Runs TBSS with specified probability model, masking specific 
#'   results and returning 1/0 if outcome of interest is/is not detected 
#' @return A binary indicator if outcome of interest is a statistical alert
#' @param df Data as formatted by ps_match() (Bernoulli) or ps_strat() (Poisson).
#' @param model Probability model, either Poisson or Bernoulli
#' @param tree_obj Hierarchical tree with columns 'parent' and 'child' as in 
#' output from tree_builder()
#' @param sig Desired level of significance.
tbss_mask <- function(df, model = c("poisson", "bernoulli"), tree_obj, sig = .05) {
  if (model == "bernoulli") {
    # format for bernoulli tbss
    treat_outcomes <- df %>%
      filter(treatment_group == 1) %>%
      count(outcome) %>%
      rename(
        "leaf" = "outcome",
        "case" = "n"
      )

    comp_outcomes <- df %>%
      filter(treatment_group == 0) %>%
      count(outcome) %>%
      rename(
        "leaf" = "outcome",
        "control" = "n"
      )

    exp_p <- df %>%
      slice(1) %>%
      pull(exp_p)

    # run tbss
    mod <- bernoulliTBSS(
      case = treat_outcomes,
      control = comp_outcomes,
      tree = tree_obj,
      p = exp_p
    )
    
  } else if (model == "poisson") {
    df_pois <- df %>% select(-ooi_label)
    mod <- poissonTBSS(
      data = df_pois,
      tree = tree_obj
    )
  }

  # determine if ooi is in result list p < .05
  ooi <- df %>%
    slice(1) %>%
    pull(ooi_label)

  if (ooi %in% (summary(mod) |> filter(pvalue < sig))$node) {
    return(1)
  } else {
    return(0)
  }
}
