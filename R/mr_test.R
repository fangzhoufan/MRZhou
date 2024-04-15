#' Return the  Cochran’s Q heterogeneity test and directional pleiotropy test
#'
#' @param dat harmonization(using the function of harmonise_data) note:the SNP>=3
#' @param full the more details of test results or none
#'
#' @return the sensitivity analysis
#' @export
#'
#' @examples
mr_test <- function(dat, full = TRUE) {
  # 检查输入数据框是否为空
  stopifnot(!is.null(dat) && nrow(dat) > 0)
  results <- data.frame()
  results2 <- data.frame()
  ## 这里我们需要dat里面最低的nsnp要大于等于3，才可以进行异质性和水平多效性检验
  # 循环处理每一对曝光和结果
  for (i in 1:length(unique(dat$exposure))) {
    strings1 <- unique(dat$exposure)[i]

    for (j in 1:length(unique(dat$outcome))) {
      strings2 <- unique(dat$outcome)[j]
      dat1 <- dat[(dat$exposure == strings1 & dat$outcome == strings2), ]

      df1 <- mr_heterogeneity(dat1)
      df1_new <- tidyr::pivot_wider(
        df1,
        names_from = "method", names_vary = "slowest",
        values_from = c("Q", "Q_df", "Q_pval")
      ) %>%
        dplyr::select(-id.exposure, -id.outcome) %>%
        as.data.frame() %>%
        dplyr::select(exposure, everything())
      df2 <- mr_pleiotropy_test(dat1)
      df2_new <- df2 %>%
        dplyr::select(-id.exposure, -id.outcome) %>%
        as.data.frame() %>%
        dplyr::select(exposure, everything()) %>%
        dplyr::rename(egger_se = se, egger_pval = pval)
      testall <- merge(df1_new, df2_new, by = c("exposure", "outcome"))

      df1_old <- df1 %>%
        dplyr::select(exposure, outcome, method, Q_pval) %>%
        as.data.frame() %>%
        dplyr::select(exposure, everything()) %>%
        dplyr::rename(pval = Q_pval)
      df2_old <- df2 %>%
        dplyr::select(exposure, outcome, pval) %>%
        as.data.frame() %>%
        dplyr::select(exposure, everything()) %>%
        mutate(method = "MR-Egger regression") %>%
        dplyr::select(exposure, outcome, method, pval)

      testallnew <- bind_rows(df1_old, df2_old)


      # 将每次迭代生成的数据框添加到结果数据框中
      results <- bind_rows(results, testall)
      results2 <- bind_rows(results2, testallnew)
    }
  }

  if (full) {
    return(results)
  } else {
    return(results2)
  }
}
