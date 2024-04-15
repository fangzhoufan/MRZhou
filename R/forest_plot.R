#' Return the forest plot
#'
#' @param data harmonization(using the function of harmonise_data)
#' @param xlim the forestploter packages param,xlim=c(0,4)
#' @param ticks_at the forestploter packages param, ticks_at=c(0,0.5,1, 2, 3)
#'
#' @return plot
#' @export
#'
#' @examples
#' library(grid)
#' library(readr)
#' library(forestploter)
forest_plot <- function(data,xlim=c(0,4),ticks_at=c(0,0.5,1, 2, 3)) {
  data<-generate_odds_ratios(data)
  data$exposure <- ifelse(data$method %in% c('Inverse variance weighted', 'Wald ratio'), data$exposure, '')
  data$outcome <- ifelse(data$method %in% c('Inverse variance weighted', 'Wald ratio'), data$outcome, '')
  data$OR <- sprintf("%.3f (%.3f to %.3f)", data$or, data$or_lci95, data$or_uci95)
  data$pval <- ifelse(data$pval < 0.001, format(data$pval, scientific = TRUE), format(round(data$pval, 3), nsmall = 3))
  data$` `<-paste(rep(" ", 22), collapse = " ")
  data<-data%>%
    dplyr::select(exposure,outcome,nsnp,method,pval,` `,OR,or,or_lci95,or_uci95)%>%
    dplyr::rename('Exposure'=exposure,
                  'Outcome'=outcome,
                  'No.of SNP'=nsnp,
                  'Method'=method,
                  'OR (95% CI)'=OR,
                  'P value'=pval
    )

  tm <- forestploter::forest_theme(
    base_size = 18,
    ci_pch = 16, ci_lty = 1, ci_lwd = 1.5, ci_col = "black", ci_Theight = 0.2,
    refline_lty = "dashed", refline_lwd = 1, refline_col = "grey20",
    xaxis_cex = 0.8,
    footnote_cex = 0.6, footnote_col = "blue"
  )

  plot <- forestploter::forest(
    data[, 1:7],
    est = data$or,
    lower = data$or_lci95,
    upper = data$or_uci95,
    ci_column = 6,
    ref_line = 1,
    xlim = xlim,
    ticks_at = ticks_at,
    footnote = "",
    theme = tm
  )

  boxcolor <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148")
  boxcolor <- boxcolor[as.numeric(as.factor(data$Method))]
  for (i in 1:nrow(data)) {
    plot <- forestploter::edit_plot(plot, col = 6, row = i, which = "ci", gp = grid::gpar(fill = boxcolor[i], fontsize = 25))
  }

  pos_bold_pval <- which(as.numeric(gsub('<', "", data$`P value`)) < 0.05)
  if (length(pos_bold_pval) > 0) {
    for (i in pos_bold_pval) {
      plot <- forestploter::edit_plot(plot, col = 5, row = i, which = "text", gp = grid::gpar(fontface = "bold"))
    }
  }

  plot <- add_border(plot, part = "header", row = 1, where = "top", gp = gpar(lwd = 2))
  plot <- add_border(plot, part = "header", row = which(data$Exposure != ''), gp = gpar(lwd = 1))
  plot <- edit_plot(plot, col = 1:ncol(data), row = 1:nrow(data), which = "text", gp = gpar(fontsize = 12))
  plot <- edit_plot(plot, col = 1:ncol(data), which = "text", hjust = unit(0.5, "npc"), part = "header",
                    x = unit(0.5, "npc"))
  plot <- edit_plot(plot, col = 1:ncol(data), which = "text", hjust = unit(0.5, "npc"),
                    x = unit(0.5, "npc"))

  return(plot)
}
