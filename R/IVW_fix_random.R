#' Return the main four MR results
#'
#' @param dat harmonization(using the function of harmonise_data)
#' @param plot the choice of showing the MR results visualization
#' @param other_method  the mr methods:"mr_egger_regression","mr_weighted_median", "mr_raps"
#' @param Random_only the choices of the only use IVW_mre
#'
#' @return
#' @export
#'
#' @examples
IVW_fix_random <- function(dat, other_method = c(
                             "mr_egger_regression",
                             "mr_weighted_median", "mr_raps"
                           ), plot = FALSE, Random_only = FALSE) {
  # 检查输入数据框是否为空
  stopifnot(!is.null(dat) && nrow(dat) > 0)
  results <- data.frame() # 初始化结果数据框
  for (i in 1:length(unique(dat$exposure))) {
    exposure_length <- length(unique(dat$exposure))
    strings1 <- unique(dat$exposure)[i]
    for (j in 1:length(unique(dat$outcome))) {
      strings2 <- unique(dat$outcome)[j]
      dat1 <- dat[(dat$exposure == strings1 & dat$outcome == strings2), ]
      outcome_length <- length(unique(dat$outcome))

      # 根据dat1的行数分类处理
      if (nrow(dat1) == 1) {
        # 对于只有一行的情况，这里填写只有一行时的处理代码
        # 例如:
        handle_single_row <- function(dat) {
          res <- mr(dat, method = c("mr_wald_ratio"))
          res <- res %>% dplyr::select(exposure, outcome, method, nsnp, b, se, pval)
          return(res)
        }
        result1 <- handle_single_row(dat1)
      } else if (nrow(dat1) == 2) {
        # 对于有两行的情况
        # 这里填写有两行时的处理代码
        # 例如:
        handle_two_rows <- function(dat, Random_only = FALSE) {
          df1 <- TwoSampleMR::mr_heterogeneity(dat)
          if (Random_only) {
            res <- mr(dat, method_list = c(
              "mr_ivw_mre",
              other_method
            ))
            res$method[res$method == "Inverse variance weighted (multiplicative random effects)"] <- "Inverse variance weighted"
          } else {
            if (df1$Q_pval < 0.05) {
              res <- mr(dat, method_list = c(
                "mr_ivw_mre",
                other_method
              ))
              res$method[res$method == "Inverse variance weighted (multiplicative random effects)"] <- "Inverse variance weighted"
            } else {
              res <- mr(dat, method_list = c(
                "mr_ivw_fe",
                other_method
              ))
              res$method[res$method == "Inverse variance weighted (fixed effects)"] <- "Inverse variance weighted"
            }
          }
          res$method[res$method == "Robust adjusted profile score (RAPS)"] <- "MR.RAPS"
          res <- res %>% dplyr::select(exposure, outcome, method, nsnp, b, se, pval)
          return(res)
        }
        result1 <- handle_two_rows(dat1)
      } else if (nrow(dat1) >= 3) {
        # 对于有三行或更多行的情况
        # 这里填写有三行或更多时的处理代码
        # 例如:
        handle_multiple_rows <- function(dat, plot = FALSE, Random_only = FALSE) {
          df1 <- TwoSampleMR::mr_heterogeneity(dat)
          if (Random_only) {
            res <- mr(dat, method_list = c(
              "mr_ivw_mre",
              other_method
            ))
            res$method[res$method == "Inverse variance weighted (multiplicative random effects)"] <- "Inverse variance weighted"
            res$method[res$method == "Robust adjusted profile score (RAPS)"] <- "MR.RAPS"

            if (plot) {
              if (res$pval[res$method == "Inverse variance weighted"] < 0.05 / (exposure_length * outcome_length)) {
                res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_egger_regression"))
                res_single$SNP[res_single$SNP == "All - Inverse variance weighted (multiplicative random effects)"] <- "All - Inverse variance weighted"
                p1 <- mr_scatter_plot(res, dat)
                a1 <- p1[[1]] + ggsci::scale_color_lancet()
                p2 <- mr_forest_plot(res_single)
                a2 <- p2[[1]] + ggsci::scale_color_lancet()
                res_loo <- mr_leaveoneout(dat, method = mr_ivw_mre)
                p3 <- mr_leaveoneout_plot(res_loo)
                a3 <- p3[[1]] + ggsci::scale_color_lancet()
                p4 <- mr_funnel_plot(res_single)
                a4 <- p4[[1]] + ggsci::scale_color_lancet()
                folder_name <- paste0(strings1, "_to_", strings2)
                dir.create(folder_name)
                ggsave(a1, file = paste0(folder_name, "/Ajusted_scatter_plot.pdf"), width = 7, height = 7)
                ggsave(a2, file = paste0(folder_name, "/Ajusted_forest_plot.pdf"), width = 7, height = 7)
                ggsave(a3, file = paste0(folder_name, "/Ajusted_leaveoneout_plot.pdf"), width = 7, height = 7)
                ggsave(a4, file = paste0(folder_name, "/Ajusted_funnel_plot.pdf"), width = 7, height = 7)
              } else if (res$pval[res$method == "Inverse variance weighted"] < 0.05) {
                res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_egger_regression"))
                res_single$SNP[res_single$SNP == "All - Inverse variance weighted (multiplicative random effects)"] <- "All - Inverse variance weighted"
                p1 <- mr_scatter_plot(res, dat)
                a1 <- p1[[1]] + ggsci::scale_color_lancet()
                p2 <- mr_forest_plot(res_single)
                a2 <- p2[[1]] + ggsci::scale_color_lancet()
                res_loo <- mr_leaveoneout(dat, method = mr_ivw_mre)
                p3 <- mr_leaveoneout_plot(res_loo)
                a3 <- p3[[1]] + ggsci::scale_color_lancet()
                p4 <- mr_funnel_plot(res_single)
                a4 <- p4[[1]] + ggsci::scale_color_lancet()
                folder_name <- paste0(strings1, "_to_", strings2)
                dir.create(folder_name)
                ggsave(a1, file = paste0(folder_name, "/scatter_plot.pdf"), width = 7, height = 7)
                ggsave(a2, file = paste0(folder_name, "/forest_plot.pdf"), width = 7, height = 7)
                ggsave(a3, file = paste0(folder_name, "/leaveoneout_plot.pdf"), width = 7, height = 7)
                ggsave(a4, file = paste0(folder_name, "/funnel_plot.pdf"), width = 7, height = 7)
              }
            }

            res <- res %>% dplyr::select(exposure, outcome, method, nsnp, b, se, pval)
          } else {
            if (df1$Q_pval[1] > 0.05 & df1$Q_pval[2] > 0.05) {
              res <- mr(dat, method_list = c(
                "mr_ivw_fe",
                other_method
              ))
              res$method[res$method == "Inverse variance weighted (fixed effects)"] <- "Inverse variance weighted"
              res$method[res$method == "Robust adjusted profile score (RAPS)"] <- "MR.RAPS"

              if (plot) {
                if (res$pval[res$method == "Inverse variance weighted"] < 0.05 / (exposure_length * outcome_length)) {
                  res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_fe", "mr_egger_regression"))
                  res_single$SNP[res_single$SNP == "All - Inverse variance weighted (fixed effects)"] <- "All - Inverse variance weighted"
                  p1 <- mr_scatter_plot(res, dat)
                  a1 <- p1[[1]] + ggsci::scale_color_lancet()
                  p2 <- mr_forest_plot(res_single)
                  a2 <- p2[[1]] + ggsci::scale_color_lancet()
                  res_loo <- mr_leaveoneout(dat, method = mr_ivw_fe)
                  p3 <- mr_leaveoneout_plot(res_loo)
                  a3 <- p3[[1]] + ggsci::scale_color_lancet()
                  p4 <- mr_funnel_plot(res_single)
                  a4 <- p4[[1]] + ggsci::scale_color_lancet()
                  folder_name <- paste0(strings1, "_to_", strings2)
                  dir.create(folder_name)
                  ggsave(a1, file = paste0(folder_name, "/Ajusted_scatter_plot.pdf"), width = 7, height = 7)
                  ggsave(a2, file = paste0(folder_name, "/Ajusted_forest_plot.pdf"), width = 7, height = 7)
                  ggsave(a3, file = paste0(folder_name, "/Ajusted_leaveoneout_plot.pdf"), width = 7, height = 7)
                  ggsave(a4, file = paste0(folder_name, "/Ajusted_funnel_plot.pdf"), width = 7, height = 7)
                } else if (res$pval[res$method == "Inverse variance weighted"] < 0.05) {
                  res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_fe", "mr_egger_regression"))
                  res_single$SNP[res_single$SNP == "All - Inverse variance weighted (fixed effects)"] <- "All - Inverse variance weighted"
                  p1 <- mr_scatter_plot(res, dat)
                  a1 <- p1[[1]] + ggsci::scale_color_lancet()
                  p2 <- mr_forest_plot(res_single)
                  a2 <- p2[[1]] + ggsci::scale_color_lancet()
                  res_loo <- mr_leaveoneout(dat, method = mr_ivw_fe)
                  p3 <- mr_leaveoneout_plot(res_loo)
                  a3 <- p3[[1]] + ggsci::scale_color_lancet()
                  p4 <- mr_funnel_plot(res_single)
                  a4 <- p4[[1]] + ggsci::scale_color_lancet()
                  folder_name <- paste0(strings1, "_to_", strings2)
                  dir.create(folder_name)
                  ggsave(a1, file = paste0(folder_name, "/scatter_plot.pdf"), width = 7, height = 7)
                  ggsave(a2, file = paste0(folder_name, "/forest_plot.pdf"), width = 7, height = 7)
                  ggsave(a3, file = paste0(folder_name, "/leaveoneout_plot.pdf"), width = 7, height = 7)
                  ggsave(a4, file = paste0(folder_name, "/funnel_plot.pdf"), width = 7, height = 7)
                }
              }
              res <- res %>% dplyr::select(exposure, outcome, method, nsnp, b, se, pval)
            } else {
              res <- mr(dat, method_list = c(
                "mr_ivw_mre",
                other_method
              ))
              res$method[res$method == "Inverse variance weighted (multiplicative random effects)"] <- "Inverse variance weighted"
              res$method[res$method == "Robust adjusted profile score (RAPS)"] <- "MR.RAPS"

              if (plot) {
                if (res$pval[res$method == "Inverse variance weighted"] < 0.05 / (exposure_length * outcome_length)) {
                  res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_egger_regression"))
                  res_single$SNP[res_single$SNP == "All - Inverse variance weighted (multiplicative random effects)"] <- "All - Inverse variance weighted"
                  p1 <- mr_scatter_plot(res, dat)
                  a1 <- p1[[1]] + ggsci::scale_color_lancet()
                  p2 <- mr_forest_plot(res_single)
                  a2 <- p2[[1]] + ggsci::scale_color_lancet()
                  res_loo <- mr_leaveoneout(dat, method = mr_ivw_mre)
                  p3 <- mr_leaveoneout_plot(res_loo)
                  a3 <- p3[[1]] + ggsci::scale_color_lancet()
                  p4 <- mr_funnel_plot(res_single)
                  a4 <- p4[[1]] + ggsci::scale_color_lancet()
                  folder_name <- paste0(strings1, "_to_", strings2)
                  dir.create(folder_name)
                  ggsave(a1, file = paste0(folder_name, "/Ajusted_scatter_plot.pdf"), width = 7, height = 7)
                  ggsave(a2, file = paste0(folder_name, "/Ajusted_forest_plot.pdf"), width = 7, height = 7)
                  ggsave(a3, file = paste0(folder_name, "/Ajusted_leaveoneout_plot.pdf"), width = 7, height = 7)
                  ggsave(a4, file = paste0(folder_name, "/Ajusted_funnel_plot.pdf"), width = 7, height = 7)
                } else if (res$pval[res$method == "Inverse variance weighted"] < 0.05) {
                  res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_egger_regression"))
                  res_single$SNP[res_single$SNP == "All - Inverse variance weighted (multiplicative random effects)"] <- "All - Inverse variance weighted"
                  p1 <- mr_scatter_plot(res, dat)
                  a1 <- p1[[1]] + ggsci::scale_color_lancet()
                  p2 <- mr_forest_plot(res_single)
                  a2 <- p2[[1]] + ggsci::scale_color_lancet()
                  res_loo <- mr_leaveoneout(dat, method = mr_ivw_mre)
                  p3 <- mr_leaveoneout_plot(res_loo)
                  a3 <- p3[[1]] + ggsci::scale_color_lancet()
                  p4 <- mr_funnel_plot(res_single)
                  a4 <- p4[[1]] + ggsci::scale_color_lancet()
                  folder_name <- paste0(strings1, "_to_", strings2)
                  dir.create(folder_name)
                  ggsave(a1, file = paste0(folder_name, "/scatter_plot.pdf"), width = 7, height = 7)
                  ggsave(a2, file = paste0(folder_name, "/forest_plot.pdf"), width = 7, height = 7)
                  ggsave(a3, file = paste0(folder_name, "/leaveoneout_plot.pdf"), width = 7, height = 7)
                  ggsave(a4, file = paste0(folder_name, "/funnel_plot.pdf"), width = 7, height = 7)
                }
              }

              res <- res %>% dplyr::select(exposure, outcome, method, nsnp, b, se, pval)
            }
          }


          return(res)
        }

        if (plot) {
          if (Random_only) {
            result1 <- handle_multiple_rows(dat1, plot = TRUE, Random_only = TRUE)
          } else {
            result1 <- handle_multiple_rows(dat1, plot = TRUE, Random_only = FALSE)
          }
        } else {
          if (Random_only) {
            result1 <- handle_multiple_rows(dat1, plot = FALSE, Random_only = TRUE)
          } else {
            result1 <- handle_multiple_rows(dat1, plot = FALSE, Random_only = FALSE)
          }
        }
      }
      # 将每次处理的结果添加到结果数据框中
      results <- rbind(results, result1)
    }
  }
  return(results)
}
