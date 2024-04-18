#' Return a circle plot
#'
#' @param a  The combine_results()  function returned a result
#' @param padj the padj plot
#'
#' @return Return a circle plot
#' @export
#'
#' @examples
#' # p+expand_limits(x = c(-100, 400), y = c(0, 4))
circle_plot <- function(a, padj = FALSE) {
  a <- a %>%
    dplyr::select(exposure, outcome, `b_Inverse variance weighted`, `pval_Inverse variance weighted`) %>%
    dplyr::rename(P = `pval_Inverse variance weighted`) %>%
    mutate(P = ifelse(`b_Inverse variance weighted` > 0, P, -P)) %>%
    dplyr::select(exposure, outcome, P)
  a <- tidyr::pivot_wider(
    a,
    names_from = "outcome", names_vary = "slowest", names_sep = ".",
    values_from = c("P")
  )

  head(a, 3)
  cirp <- melt(a)
  head(cirp, 3)
  res <- data.frame(exposure = unique(a$exposure))
  res$ang <- seq(
    from = (360 / nrow(res)) / 1.5,
    to = (1.5 * (360 / nrow(res))) - 360,
    length.out = nrow(res)
  ) + 80
  res$hjust <- 0
  res$hjust[which(res$ang < -90)] <- 1
  res$ang[which(res$ang < -90)] <- (180 + res$ang)[which(res$ang < -90)]
  cirp$var <- rep(1:length(unique(a$exposure)), length(colnames(a)) - 1)

  if (padj) {
    bon <- (0.05 / nrow(a))
    cuts <- c(-1, -0.05, -bon, 0, bon, 0.05, 1)
    bon <- ifelse(bon < 0.01, format(bon, scientific = TRUE, digits = 2), format(floor(bon * 100) / 100, nsmall = 2))
    labels <- c("NC", "p<0.05(Low risk)", paste0("p<", bon, "(Low risk)"), paste0("p<", bon, "(High risk)"), "p<0.05(High risk)", "NC")
    colors <- c("red1", "orangered", "skyblue1", "royalblue1", "blue4")
    cirp$value_group <- cut(cirp$value, breaks = cuts, labels = labels)
    cirp$value_group <- factor(cirp$value_group, levels = c(paste0("p<", bon, "(High risk)"), "p<0.05(High risk)", "NC", "p<0.05(Low risk)", paste0("p<", bon, "(Low risk)")))
  } else {
    cuts <- c(-1, -0.05, 0, 0.05, 1)
    labels <- c("NC", "p<0.05(Low risk)", "p<0.05(High risk)", "NC")
    colors <- c("red1", "skyblue1", "blue4")
    cirp$value_group <- cut(cirp$value, breaks = cuts, labels = labels)
    cirp$value_group <- factor(cirp$value_group, levels = c("p<0.05(High risk)", "NC", "p<0.05(Low risk)"))
  }
  strings <- colnames(a)[-1]
  length(strings)
  num_strings <- length(strings)
  p <- ggplot(cirp, aes(x = var, y = variable)) +
    geom_tile(aes(fill = value_group), color = "white") +
    scale_fill_manual(values = colors) +
    scale_y_discrete(expand = expansion(mult = c(2, 0))) +
    coord_polar(theta = "x") +
    theme_void() +
    geom_text(
      data = res,
      aes(
        x = as.numeric(rownames(res)),
        y = num_strings + 0.6,
        label = exposure, angle = ang, hjust = hjust
      ),
      size = 3
    ) +
    xlim(0, nrow(res) * 8 / 7) +
    theme(legend.position = c(0.5, 0.5)) +
    guides(fill = guide_legend(title = "Risks(IVW method)")) +
    expand_limits(x = c(-100, 200), y = c(0, 5))
  p <- p + lapply(1:num_strings, function(i) {
    geom_text(
      data = cirp[which(cirp$variable == strings[i]), ],
      aes(x = 0.5, y = i, label = variable),
      position = position_dodge(width = 0.5),
      hjust = 1, vjust = 0.5, size = 3
    )
  })
  p
  return(p)
}
