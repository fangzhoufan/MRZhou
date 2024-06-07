#' SMR data prepration
#'
#' @param gwas the read_vcf results dataframe
#' @param snp_col SNP
#' @param A1_col effect_allele
#' @param A2_col other_allele
#' @param freq_col eaf
#' @param b_col beta
#' @param se_col se
#' @param p_col pval
#' @param N_col samplesize
#' @param filename filename,eg:'LUAD.ma'
#'
#' @return
#' @export
#'
#' @examples
smr_data <- function(gwas,
                     snp_col = "SNP",
                     A1_col = "effect_allele",
                     A2_col = "other_allele",
                     freq_col = "eaf",
                     b_col = "beta",
                     se_col = "se",
                     p_col = "pval",
                     N_col = "samplesize",
                     filename) {
  gwas <- gwas %>%
    dplyr::rename(
      `SNP` = snp_col,
      `A1` = A1_col,
      `A2` = A2_col,
      `freq` = freq_col,
      `b` = b_col,
      `se` = se_col,
      `p` = p_col,
      `N` = N_col
    ) %>%
    dplyr::filter(str_detect(SNP, "^rs")) %>%
    dplyr::select(SNP, A1, A2, freq, b, se, p, N) %>%
    dplyr::filter(!str_detect(SNP, ",")) %>%
    dplyr::distinct(SNP)

  write.table(gwas, paste0(filename), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  return(gwas)
}
