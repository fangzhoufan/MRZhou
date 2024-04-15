#' Return exposure when using the finngen database
#'
#' @param trait the name of exposure, eg: osteonecrosis
#' @param filepath the exposure original file path
#' @param samplesize  the exposure samplesize
#' @param p_threshold  the p_threshold = 5e-8
#'
#' @return a dataframe of finngen exposure IVs
#' @export
#'
#' @examples
read_finngen_exposure <- function(trait, filepath, samplesize, p_threshold = 5e-8) {
  ## 读取数据
  expo_rt <- read.csv(filepath, header = TRUE)

  ## 添加样本大小和trait列
  expo_rt$n <- samplesize
  expo_rt$trait <- trait

  ## 过滤p值小于指定阈值的行
  expo_rt <- expo_rt[expo_rt$pval < p_threshold, ]

  ## 检查是否有数据满足条件
  if (nrow(expo_rt) == 0) {
    message("No data with p-value <", p_threshold, " found.")
    return(data.frame())
  }

  ## 写入文件
  write.csv(expo_rt, file = "d.csv")

  ## 读取处理后的文件
  df <- read_exposure_data(
    filename = "d.csv",
    sep = ",",
    samplesize_col = "n",
    phenotype_col = "trait",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval",
    gene_col = "nearest_genes",
    chr_col = "#chrom",
    pos_col = "pos"
  )
  df <- clump_data(df, clump_r2 = 0.001, clump_kb = 10000)

  ## 如果没有数据，提示并返回空的数据框
  if (nrow(df) == 0) {
    message("No data found after processing.")
    return(data.frame())
  }
  return(df)
}
