#' Return the MR-Presso results of dat (after harmonization)
#'
#' @param dat harmonization(using the function of harmonise_data)
#'
#' @return the value of MR-Presso results
#'        (exposure,outcome,Presso.Pvalue,Presso.RSSobs,b.MR.presso,se.MR.presso,pval.MR.presso)
#'
#' @export
#'
#' @examples
MRpresso <- function(dat) {
  # 创建一个空的结果数据框
  results <- data.frame(Exposure = character(),
                        Outcome = character(),
                        MR.presso.global.RSSobs = character(),
                        MR.presso.global.pvalue = character(),
                        b.MR.presso = numeric(),
                        se.MR.presso = numeric(),
                        pval.MR.presso = numeric(),
                        stringsAsFactors = FALSE)

  # 循环处理每一对曝露和结果
  for (i in 1:length(unique(dat$exposure))) {
    strings1 <- unique(dat$exposure)[i]

    for (j in 1:length(unique(dat$outcome))) {
      strings2 <- unique(dat$outcome)[j]
      dat1 <- dat[(dat$exposure == strings1 & dat$outcome == strings2), ]

      # 尝试运行 MR-PRESSO 分析
      res_presso <- tryCatch({
        TwoSampleMR::run_mr_presso(dat1)
      }, error = function(e) {
        print(paste0(strings1, ' to ', strings2, ' Not found'))
        NA
      })

      # 处理 MR-PRESSO 结果
      if (!is.na(res_presso)) {
        Presso.Pvalue <- as.character(res_presso[[1]][["MR-PRESSO results"]][["Global Test"]][["Pvalue"]])
        Presso.RSSobs <- as.character(res_presso[[1]][["MR-PRESSO results"]][["Global Test"]][["RSSobs"]])
        b.MR.presso <- res_presso[[1]][["Main MR results"]][["Causal Estimate"]][1]
        se.MR.presso <- res_presso[[1]][["Main MR results"]][["Sd"]][1]
        pval.MR.presso <- res_presso[[1]][["Main MR results"]][["P-value"]][1]
      } else {
        # 如果发生错误，将所有变量设置为NA或其他值
        Presso.Pvalue <- NA
        Presso.RSSobs <- NA
        b.MR.presso <- NA
        se.MR.presso <- NA
        pval.MR.presso <- NA
      }

      # 创建包含结果的数据框
      result <- data.frame(Exposure = strings1,
                           Outcome = strings2,
                           MR.presso.global.RSSobs = Presso.RSSobs,
                           MR.presso.global.pvalue = Presso.Pvalue,
                           b.MR.presso = b.MR.presso,
                           se.MR.presso = se.MR.presso,
                           pval.MR.presso = pval.MR.presso,
                           stringsAsFactors = FALSE)

      # 将结果追加到结果数据框中
      results <- rbind(results, result)
    }
  }

  # 返回结果数据框
  return(results)
}
