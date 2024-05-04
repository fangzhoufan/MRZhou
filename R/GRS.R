#' Return the genetic risk scores
#'
#' @param dat the harmonization data
#'
#' @return GRS dataframe
#' @export
#'
#' @examples
GRS <- function(dat) {
  exposures <- unique(dat$exposure)
  outcomes <- unique(dat$outcome)
  grs_data_all <- data.frame()  # 创建一个空数据框来存储所有结果
  for (exp in exposures) {
    for (out in outcomes) {
      grs_dat <- dat %>%
        dplyr::filter(exposure == exp) %>%
        dplyr::filter(outcome == out) %>%
        dplyr::filter(mr_keep == 'TRUE')
      if (nrow(grs_dat) > 0) {  # 确保筛选后的数据框不为空
        grs <- grs.summary((grs_dat$beta.exposure / grs_dat$se.exposure),
                           grs_dat$beta.outcome, grs_dat$se.outcome, n = nrow(grs_dat))
        OR <- round(exp(grs$ahat), 6)
        CI_lower <- round(exp(grs$ahat - 1.96 * grs$aSE), 6)
        CI_upper <- round(exp(grs$ahat + 1.96 * grs$aSE), 6)
        ORs <- paste0(OR, ' (', CI_lower, '-', CI_upper, ')')
        p_val <- grs$pval
        if (p_val < 0.01) {
          p_val_formatted <- sprintf("%.2e", p_val)
        } else {
          p_val_formatted <- round(p_val, 2)
        }
        grs_data <- data.frame(Exposure = exp, Outcome = out, OR = ORs, p = p_val_formatted)
        grs_data <- grs_data %>% dplyr::rename('OR(95%CI)' = OR, 'P value' = p)
        grs_data_all <- rbind(grs_data_all, grs_data)  # 将每次计算的结果添加到总数据框中
      }
    }
  }

  return(grs_data_all)
}
