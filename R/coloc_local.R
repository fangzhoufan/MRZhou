#' Prepare datasets to Bayesian colocalisation analysis using Bayes Factors
#'
#' @param eqtl the dataframe of eqtl
#' @param gwas the dataframe of gwas
#' @param beta_col_eqtl regression coefficient for each SNP from dataset1
#' @param se_col_eqtl variance of beta from dataset1
#' @param snp_col_eqtl a character vector of snp ids, optional. It will be used to merge dataset1 and dataset2 and will be retained in the results.
#' @param pval_col_eqtl P-values for each SNP in dataset1
#' @param MAF_col_eqtl minor allele frequency of the variants
#' @param N_eqtl Number of samples in dataset 1
#' @param type_eqtl the type of data in dataset 1 - either "quant" or "cc" to denote quantitative or case-control
#' @param s_eqtl for a case control dataset, the proportion of samples in dataset 1 that are cases
#' @param sdY_eqtl for a quantitative trait, the population standard deviation of the trait. if not given, it can be estimated from the vectors of varbeta and MAF
#' @param beta_col_gwas regression coefficient for each SNP from gwas
#' @param se_col_gwas variance of beta from gwas
#' @param snp_col_gwas a character vector of snp ids, optional. It will be used to merge dataset1 and dataset2 and will be retained in the results.
#' @param pval_col_gwas P-values for each SNP in dataset1
#' @param MAF_col_gwas minor allele frequency of the variants
#' @param N_gwas Number of samples in dataset 2
#' @param type_gwas the type of data in dataset 2 - either "quant" or "cc" to denote quantitative or case-control
#' @param s_gwas for a case control dataset, the proportion of samples in dataset 2 that are cases
#' @param sdY_gwas for a quantitative trait, the population standard deviation of the trait. if not given, it can be estimated from the vectors of varbeta and MAF
#'
#' @return the list of dataset 1 and 2
#' @export
#'
#' @examples
#' # also can:
#' # check_dataset(result$dataset1)
#' # check_dataset(result$dataset2)
#' # my.res <- coloc.abf(result$dataset1, result$dataset2)
#'
#' # only using the pval and SNP and maf
#' # result_coloc <- coloc.abf(dataset1=list(snp=result$dataset1$snp,pvalues=result$dataset1$pvalues,type="cc", s=0.001, N=359399), dataset2=list(snp=result$dataset2$snp,pvalues=result$dataset2$pvalues, type="quant", N=31684), MAF=result$dataset1$MAF)
coloc_local <- function(eqtl, gwas,
                        beta_col_eqtl = "beta",
                        se_col_eqtl = "varbeta",
                        snp_col_eqtl = "snp",
                        pval_col_eqtl = "pval",
                        MAF_col_eqtl = "MAF",
                        N_eqtl = NULL,
                        type_eqtl = NULL,
                        s_eqtl = NULL,
                        sdY_eqtl = NULL,
                        beta_col_gwas = "beta",
                        se_col_gwas = "varbeta",
                        snp_col_gwas = "snp",
                        pval_col_gwas = "pval",
                        MAF_col_gwas = "MAF",
                        N_gwas = NULL,
                        type_gwas = NULL,
                        s_gwas = NULL,
                        sdY_gwas = NULL) {
  # 检查列名是否存在
  check_columns_exist <- function(df, col_names) {
    if (any(!col_names %in% names(df))) {
      missing <- col_names[!col_names %in% names(df)]
      stop(paste("The following columns are missing in the data frame:", paste(missing, collapse = ", ")))
    }
  }

  # 检查eqtl数据框的列名
  check_columns_exist(eqtl, c(beta_col_eqtl, snp_col_eqtl, pval_col_eqtl, MAF_col_eqtl))
  # 检查gwas数据框的列名
  check_columns_exist(gwas, c(beta_col_gwas, snp_col_gwas, pval_col_gwas, MAF_col_gwas))

  # 重命名列
  eqtl_renamed <- eqtl %>%
    dplyr::rename(
      `beta` = beta_col_eqtl,
      `varbeta` = se_col_eqtl,
      `snp` = snp_col_eqtl,
      `pval` = pval_col_eqtl,
      `MAF` = MAF_col_eqtl
    )

  gwas_renamed <- gwas %>%
    dplyr::rename(
      `beta` = beta_col_gwas,
      `varbeta` = se_col_gwas,
      `snp` = snp_col_gwas,
      `pval` = pval_col_gwas,
      `MAF` = MAF_col_gwas
    )

  input <- merge(eqtl_renamed, gwas_renamed, by = "snp", all = FALSE, suffixes = c("_eqtl", "_gwas"))
  input <- input[!duplicated(input$snp) & input$MAF_eqtl != 0 & input$MAF_gwas != 0, ]

  # 构建dataset1列表，仅包含非NULL值
  dataset1 <- list(
    snp = input$snp,
    beta = input$beta_eqtl,
    varbeta = input$varbeta_eqtl,
    pvalues = input$pval_eqtl,
    MAF = input$MAF_eqtl,
    type = if (!is.null(type_eqtl)) type_eqtl else NULL,
    s = if (!is.null(s_eqtl)) s_eqtl else NULL,
    N = if (!is.null(N_eqtl)) N_eqtl else NULL,
    sdY = if (!is.null(sdY_eqtl)) sdY_eqtl else NULL
  )

  # 移除dataset1中的NULL条目
  dataset1 <- dataset1[!sapply(dataset1, is.null)]

  # 构建dataset2列表，类似dataset1的处理
  # ... [类似dataset1的处理]
  dataset2 <- list(
    snp = input$snp,
    beta = input$beta_gwas,
    varbeta = input$varbeta_gwas,
    pvalues = input$pval_gwas,
    MAF = input$MAF_gwas,
    type = if (!is.null(type_gwas)) type_gwas else NULL,
    s = if (!is.null(s_gwas)) s_gwas else NULL,
    N = if (!is.null(N_gwas)) N_gwas else NULL,
    sdY = if (!is.null(sdY_gwas)) sdY_gwas else NULL
  )
  dataset2 <- dataset2[!sapply(dataset2, is.null)]

  dataset <- list(dataset1 = dataset1, dataset2 = dataset2)
  dataset <- compact(dataset)
  # 返回合并后的数据框
  return(dataset)
}
