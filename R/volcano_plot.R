#' Return a volcano plot
#'
#' @param mydata The combine_results()  function returned a result
#' @param adjust The choose of 'BH','FDR','None'
#'
#' @return Return a volcano plot
#' @export
#'
#' @examples
#' library(ggrepel)
#' library(gt)
volcano_plot <- function(mydata,adjust='BH') {
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
  if(adjust=='BH'){
    mydata1 <- mydata %>%
      mutate(expression = case_when(
        OR >= 1 & Pval < 0.05 ~ "High risk",
        OR <= 1 & Pval < 0.05 ~ "Low risk",
        TRUE ~ "Insignificance"
      ))
    # 选择前十个数据
    top_10 <- bind_rows(
      mydata1 %>%
        filter(Pval < 0.05) %>%
        arrange(Pval, desc(abs(OR))) %>%
        head()
    )
    adjustedp <- 0.05 / nrow(mydata)
    boni <- mydata %>%
      filter(Pval < adjustedp)
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
    volc_plot2 <- volc_plot1 +
      geom_label_repel(data = top_10,aes(OR, -log10(Pval), label = ID), size = 2)
  }else if(adjust=='FDR'){
    mydata$padj <- p.adjust(mydata$Pval,method='fdr')
    mydata1 <- mydata %>%
      mutate(expression = case_when(
        OR >= 1 & padj < 0.05 ~ "High risk",
        OR <= 1 & padj < 0.05 ~ "Low risk",
        TRUE ~ "Insignificance"
      ))
    # 选择前十个数据
    top_10 <- bind_rows(
      mydata1 %>%
        filter(padj < 0.05) %>%
        arrange(padj, desc(abs(OR))) %>%
        head()
    )
    volc_plot1 <- ggplot(mydata1, aes(OR, -log10(padj))) +
      geom_point(aes(color = expression), size = 1.4) +
      xlab("OR") +
      ylab("-log10(Padj)") +
      scale_x_continuous(limits = c(min(mydata$OR) - 0.1, max(mydata$OR) + 0.1)) +
      scale_color_manual(values = c("red1", "grey", "skyblue")) +
      guides(color = guide_legend(title = "Risks(IVW method)")) +
      # 添加分界线
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
    volc_plot2 <- volc_plot1 +
      geom_label_repel(data = top_10,aes(OR, -log10(padj), label = ID), size = 2)
  }else if(adjust=='None'){
    mydata1 <- mydata %>%
      mutate(expression = case_when(
        OR >= 1 & Pval < 0.05 ~ "High risk",
        OR <= 1 & Pval < 0.05 ~ "Low risk",
        TRUE ~ "Insignificance"
      ))
    # 选择前十个数据
    top_10 <- bind_rows(
      mydata1 %>%
        filter(Pval < 0.05) %>%
        arrange(Pval, desc(abs(OR))) %>%
        head()
    )
    volc_plot1 <- ggplot(mydata1, aes(OR, -log10(Pval))) +
      geom_point(aes(color = expression), size = 1.4) +
      xlab("OR") +
      ylab("-log10(P-value)") +
      scale_x_continuous(limits = c(min(mydata$OR) - 0.1, max(mydata$OR) + 0.1)) +
      scale_color_manual(values = c("red1", "grey", "skyblue")) +
      guides(color = guide_legend(title = "Risks(IVW method)")) +
      # 添加分界线
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
    volc_plot2 <- volc_plot1 +
      geom_label_repel(data = top_10,aes(OR, -log10(Pval), label = ID), size = 2)
  }

  return(volc_plot2)
}
