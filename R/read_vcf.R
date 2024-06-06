#' Read IEU GWAS database vcf files (original)
#'
#' @param path The path of vcf files
#'
#' @return  a dataframe
#' @export
#'
#' @examples
#' # library(gwasvcf)
#' # library(gwasglue)
#' # library(VariantAnnotation)
#' # library(TwoSampleMR)
read_vcf <- function(path) {
  vcf_exposure <- VariantAnnotation::readVcf(path)
  # exposure_pval_filter <- query_gwas(vcf = vcf_exposure, pval = 1)
  exposure_pval_filter <- gwasvcf_to_TwoSampleMR(vcf = vcf_exposure)
  colnames(exposure_pval_filter) <- gsub("\\.exposure", "", colnames(exposure_pval_filter))
  gc()
  return(exposure_pval_filter)
}
