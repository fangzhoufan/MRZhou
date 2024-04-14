#' Return the IVs (F>10)
#'
#' @param dat the exposure or harmonization between exposure and outcome
#' @param F_value the F statistic = 10
#'
#' @return Calculate the F statistic for the exposure IVs (F>10)
#' @export
#'
#' @examples
get_f <- function(dat, F_value = 10) {
  log <- is.na(dat$eaf.exposure)
  log <- unique(log)
  if (length(log) == 1) {
    if (log == TRUE) {
      print("数据不包含eaf，无法计算F统计量")
      return(dat)
    }
  }
  if (is.null(dat$beta.exposure[1]) == T || is.na(dat$beta.exposure[1]) == T) {
    print("数据不包含beta，无法计算F统计量")
    return(dat)
  }
  if (is.null(dat$se.exposure[1]) == T || is.na(dat$se.exposure[1]) == T) {
    print("数据不包含se，无法计算F统计量")
    return(dat)
  }
  if (is.null(dat$samplesize.exposure[1]) == T || is.na(dat$samplesize.exposure[1]) == T) {
    print("数据不包含samplesize(样本量)，无法计算F统计量")
    return(dat)
  }

  if ("FALSE" %in% log && is.null(dat$beta.exposure[1]) == F && is.na(dat$beta.exposure[1]) == F && is.null(dat$se.exposure[1]) == F && is.na(dat$se.exposure[1]) == F && is.null(dat$samplesize.exposure[1]) == F && is.na(dat$samplesize.exposure[1]) == F) {
    R2 <- (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$beta.exposure^2)) / ((2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$beta.exposure^2)) + (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * (dat$se.exposure^2) * dat$samplesize.exposure))
    F <- (dat$samplesize.exposure - 2) * R2 / (1 - R2)
    dat$R2 <- R2
    dat$F <- F
    dat <- subset(dat, F > F_value)
    return(dat)
  }
}
