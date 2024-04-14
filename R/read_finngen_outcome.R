#' Return outcome when using the finngen database
#'
#' @param exp the exposure IVs selection
#' @param trait the name of outcome,eg:osteonecrosis
#' @param filepath the outcome original file path
#' @param samplesize the outcome samplesize
#'
#' @return finngen outcome
#' @export
#'
#' @examples
read_finngen_outcome <- function(exp, trait, filepath, samplesize) {
  out <- fread(filepath, header = T, data.table = F)
  out$trait <- trait
  out$n <- samplesize
  out <- dplyr::rename(out, SNP = rsids)
  out <- merge(exp, out, by.x = "SNP", by.y = "SNP")
  write.csv(out, file = "d.csv")
  out <- TwoSampleMR::read_outcome_data(
    snps = exp$SNP,
    filename = "d.csv",
    sep = ",",
    samplesize_col = "n",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval", gene_col = "nearest_genes",
    phenotype_col = "trait",
    chr_col = "#chrom",
    pos_col = "pos"
  )
  return(out)
}

library(data.table)
