#' Combination the result,test,presso results
#'
#' @param result The IVW_fix_random() function returned a result
#' @param test The mr_test()  function returned a result
#' @param presso The MRpresso()  function returned a result
#'
#' @return Combination the result,test,presso results
#' @export
#'
#' @examples
combine_results <- function(result = NULL, test = NULL, presso = NULL) {
  if (!is.null(result) & !is.null(test) & !is.null(presso)) {
    res <- tidyr::pivot_wider(
      result,
      names_from = "method", names_sep = "_",
      values_from = c("b", "se", "pval")
    )
    res_all <- merge(res, presso, by = c("exposure", "outcome"), all = TRUE)
    res_all <- merge(res_all, test, by = c("exposure", "outcome"), all = TRUE)
  } else if (!is.null(result) & !is.null(presso)) {
    res <- tidyr::pivot_wider(
      result,
      names_from = "method", names_sep = "_",
      values_from = c("b", "se", "pval")
    )
    res_all <- merge(res, presso, by = c("exposure", "outcome"), all = TRUE)
    res_all <- res_all %>%
      dplyr::select(-c("MR.presso.global.RSSobs", "MR.presso.global.pvalue"))
    res_all <- tidyr::pivot_longer(res_all,
      cols = starts_with("b_") | starts_with("se_") | starts_with("pval_"),
      names_to = c(".value", "method"),
      names_sep = "_"
    ) %>%
      na.omit()
  } else if (!is.null(result) & !is.null(test)) {
    res <- tidyr::pivot_wider(
      result,
      names_from = "method", names_sep = "_",
      values_from = c("b", "se", "pval")
    )
    res_all <- merge(res, test, by = c("exposure", "outcome"), all = TRUE)
  } else if (!is.null(test) & !is.null(presso)) {
    res_all <- merge(test, presso, by = c("exposure", "outcome"), all = TRUE)
  }

  return(res_all)
}
