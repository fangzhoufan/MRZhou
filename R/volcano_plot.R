#' Return a volcano plot
#'
#' @param mydata The combine_results()  function returned a result
#'
#' @return Return a volcano plot
#' @export
#'
#' @examples
#' library(ggrepel)
#' library(gt)
volcano_plot <- function(mydata) {
  # 数据处理
  mydata <- mydata %>%
    dplyr::select(exposure, `b_Inverse variance weighted`, `pval_Inverse variance weighted`) %>%
    dplyr::mutate(OR = exp(`b_Inverse variance weighted`)) %>%
    dplyr::rename(ID = exposure, Pval = `pval_Inverse variance weighted`) %>%
    dplyr::select(ID, OR, Pval)
  # 创建火山图
  volc_plot <- ggplot(mydata, aes(OR, -log10(Pval))) +
    geom_point(size = 0.4) +
    xlab("OR") +
    ylab("-log10(P-value)") +
    scale_x_continuous(limits = c(min(mydata$OR) - 0.1, max(mydata$OR) + 0.1)) +
    theme_minimal()
  # 数据处理
  adjustedp <- 0.05 / nrow(mydata)
  boni <- mydata %>%
    filter(Pval < adjustedp)

  mydata1 <- mydata %>%
    mutate(expression = case_when(
      OR >= 1 & Pval < 0.05 ~ "High risk",
      OR <= 1 & Pval < 0.05 ~ "Low risk",
      TRUE ~ "Insignificance"
    ))

  # 创建着色的火山图
  volc_plot1 <- ggplot(mydata1, aes(OR, -log10(Pval))) +
    geom_point(aes(color = expression), size = 1.4) +
    xlab("OR") +
    ylab("-log10(P-value)") +
    scale_x_continuous(limits = c(min(mydata$OR) - 0.1, max(mydata$OR) + 0.1)) +
    scale_color_manual(values = c("red1", "grey", "skyblue")) +
    guides(color = guide_legend(title = "Risks(IVW method)")) +
    # 添加分界线
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(adjustedp), linetype = "dashed", color = "black")


  # 选择前十个数据
  top_10 <- bind_rows(
    mydata1 %>%
      filter(Pval < 0.05) %>%
      arrange(Pval, desc(abs(OR))) %>%
      head()
  )

  # 在火山图上添加标签
  volc_plot2 <- volc_plot1 +
    geom_label_repel(
      data = top_10,
      aes(OR, -log10(Pval), label = ID),
      size = 2
    )

  plot_list <- list(
    "volcano_plot" = volc_plot,
    "colored_volcano_plot" = volc_plot2,
    "head_data" = head(mydata) %>% gt(),
    "head_processed_data" = head(mydata1) %>% gt(),
    "top_10_data" = top_10 %>% gt()
  )

  return(plot_list[["colored_volcano_plot"]])
}
