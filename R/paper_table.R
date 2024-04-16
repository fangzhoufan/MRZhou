#' Return the main MR analysis results including:TableS1-3
#'
#' @param exp Exposure instrumental variable
#' @param out Outcome
#' @param dat harmonization(using the function of harmonise_data)
#' @param result The IVW_fix_random() function returned a result
#' @param test The mr_test()  function returned a result
#' @param presso The MRpresso()  function returned a result
#'
#' @return a paper needing marital
#' @export
#'
#' @examples
paper_table <- function(exp, out, dat, result, test, presso) {
  TableS1 <- exp %>%
    dplyr::select(
      exposure, SNP, chr.exposure, pos.exposure, effect_allele.exposure, other_allele.exposure,
      eaf.exposure, beta.exposure, se.exposure, pval.exposure, R2, F
    ) %>%
    dplyr::rename(
      Exposure = exposure,
      SNP = SNP,
      Chr = chr.exposure,
      Pos = pos.exposure,
      Effect = effect_allele.exposure,
      Allele = other_allele.exposure,
      "Effect Allele Frequency" = eaf.exposure,
      Beta = beta.exposure,
      SE = se.exposure,
      "P value" = pval.exposure
    )
  TableS2 <- dat
  TableS3 <- combin_results(result, test, presso)
  list <- list(
    TableS1 = TableS1,
    TableS2 = TableS2,
    TableS3 = TableS3
  )
  return(list)
}
