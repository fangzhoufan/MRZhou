#' Generates a data frame of pairwise linkage disequilibrium statistics
#'
#' @param dat a dataframe  including SNP column such as:variants while all SNP must be on the same chromosome, using an rsID or chromosome coordinate
#' @param genome_build Choose between one of the three options...'grch37' for genome build GRCh37 (hg19), 'grch38' for GRCh38 (hg38), or 'grch38_high_coverage' for GRCh38 High Coverage (hg38) 1000 Genome Project data sets. Default is GRCh37 (hg19).
#'
#' @return a data frame
#' @export
#'
#' @examples
calculate_LD <- function(dat, genome_build='grch37') {
  unique_exposures <- unique(dat$exposure)
  result_list <- list()
  for (exposure in unique_exposures) {
    snps_subset <- dat$SNP[dat$exposure == exposure]
    result <- LDmatrix(snps = snps_subset,
                       r2d = "r2",
                       token = '57275ed3e372',
                       genome_build = genome_build)
    result_df <- as.data.frame(result)
    colnames(result_df)[1] <- exposure
    result_list[[exposure]] <- result_df
  }
  return(result_list)
}
