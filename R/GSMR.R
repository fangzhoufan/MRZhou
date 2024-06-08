#' Return the GSMR analysis
#'
#' @param dat the harmonization data
#' @param plot the choice of showing the picture
#' @param folder the folder names
#' @param LDmatrix  The choice of LD correlation matrix
#'
#' @return Return the GSMR analysis dataframe
#' @export
#'
#' @examples # https://yanglab.westlake.edu.cn/software/gsmr/
GSMR <- function(dat, plot = FALSE, folder = NULL, LDmatrix = FALSE) {
  if (!is.null(folder)) {
    dir.create(folder, showWarnings = FALSE)
  }
  exposures <- unique(dat$exposure)
  outcomes <- unique(dat$outcome)
  gsmr_data_all <- data.frame()
  for (exp in exposures) {
    for (out in outcomes) {
      gsmr_data <- dat %>%
        dplyr::filter(exposure == exp) %>%
        dplyr::filter(outcome == out) %>%
        dplyr::filter(mr_keep == "TRUE")
      if (nrow(gsmr_data) > 0) { # 确保筛选后的数据框不为空
        gsmr_data <- gsmr_data %>%
          dplyr::select(
            "SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure",
            "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
            "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome"
          ) %>%
          dplyr::rename(
            "SNP" = SNP,
            "a1" = effect_allele.exposure,
            "a2" = other_allele.exposure,
            "a1_freq" = eaf.exposure,
            "bzx" = beta.exposure,
            "bzx_se" = se.exposure,
            "bzx_pval" = pval.exposure,
            "bzx_n" = samplesize.exposure,
            "bzy" = beta.outcome,
            "bzy_se" = se.outcome,
            "bzy_pval" = pval.outcome,
            "bzy_n" = samplesize.outcome
          )
        gsmr_data2 <- gsmr_data
        perform_gsmr_analysis <- function(gsmr_data, ldrho, snp_coeff_id) {
          bzx <- gsmr_data$bzx # SNP effects on the risk factor
          bzx_se <- gsmr_data$bzx_se # standard errors of bzx
          bzx_pval <- gsmr_data$bzx_pval # p-values for bzx
          bzy <- gsmr_data$bzy # SNP effects on the disease
          bzy_se <- gsmr_data$bzy_se # standard errors of bzy
          bzy_pval <- gsmr_data$bzy_pval # p-values for bzy
          n_ref <- 503 # Sample size of the reference sample used to calculate the LD matrix
          gwas_thresh <- 5e-8 # GWAS threshold to select SNPs as the instruments for the GSMR analysis
          multi_snps_heidi_thresh <- 0.01 # p-value threshold for multi-SNP-based global HEIDI-outlier analysis
          nsnps_thresh <- 2 # the minimum number of instruments required for the GSMR analysis
          heidi_outlier_flag <- T # flag for HEIDI-outlier analysis
          ld_r2_thresh <- 0.05 # LD r2 threshold to remove SNPs in high LD
          ld_fdr_thresh <- 0.05 # FDR threshold to remove the chance correlations between the SNP instruments
          gsmr2_beta <- 1 # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development
          gsmr_results <- gsmr2::gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta) # GSMR analysis
          return(gsmr_results)
        }
        if (LDmatrix) {
          ldrho <- ieugwasr::ld_matrix(gsmr_data$SNP, with_alleles = FALSE)
          # Estimate LD correlation matrix using R
          snp_coeff_id <- colnames(ldrho)
          snp_coeff <- ldrho

          # Match the SNP genotype data with the summary data
          snp_id <- Reduce(intersect, list(gsmr_data$SNP, snp_coeff_id))
          gsmr_data <- gsmr_data[match(snp_id, gsmr_data$SNP), ]
          snp_order <- match(snp_id, snp_coeff_id)
          snp_coeff_id <- snp_coeff_id[snp_order]
          snp_coeff <- snp_coeff[, snp_order]

          # Calculate the LD correlation matrix
          ldrho <- cor(snp_coeff)

          # Check the size of the correlation matrix and double-check if the order of the SNPs in the LD correlation matrix is consistent with that in the GWAS summary data
          colnames(ldrho) <- rownames(ldrho) <- snp_coeff_id
          gsmr_results <- tryCatch(
            {
              perform_gsmr_analysis(gsmr_data, ldrho, snp_coeff_id)
            },
            error = function(e) {
              gsmr_data <- gsmr_data2
              custom_error_message <- "LD矩阵时发生错误。我们进行假设LD没有相关性"
              cat("计算", exp, " ~ ", out, "：", custom_error_message, "\n")
              cat("原始错误信息：", e$message, "\n")

              # 清除第一部分运行时创建的变量
              rm(list = c("ldrho", "snp_coeff_id", "snp_coeff", "snp_id", "snp_order"))
              ldrho <- diag(nrow(gsmr_data))
              colnames(ldrho) <- rownames(ldrho) <- snp_coeff_id <- snpid <- as.character(gsmr_data$SNP)
              perform_gsmr_analysis(gsmr_data, ldrho, snp_coeff_id)
            }
          )
        } else {
          ldrho <- diag(nrow(gsmr_data))
          colnames(ldrho) <- rownames(ldrho) <- snp_coeff_id <- snpid <- as.character(gsmr_data$SNP)
          gsmr_results <- perform_gsmr_analysis(gsmr_data, ldrho, snp_coeff_id)
        }



        filtered_index <- gsmr_results$used_index

        b <- round(gsmr_results$bxy, 6)
        se <- round(gsmr_results$bxy_se, 6)
        OR <- round(exp(gsmr_results$bxy), 6)
        CI_lower <- round(exp(gsmr_results$bxy - 1.96 * gsmr_results$bxy_se), 6)
        CI_upper <- round(exp(gsmr_results$bxy + 1.96 * gsmr_results$bxy_se), 6)
        ORs <- paste0(OR, " (", CI_lower, "-", CI_upper, ")")
        p_val <- gsmr_results$bxy_pval
        if (p_val < 0.01) {
          p_val_formatted <- sprintf("%.2e", p_val)
        } else {
          p_val_formatted <- round(p_val, 2)
        }
        if (plot && !is.null(folder)) {
          effect_col <- colors()[75]
          data <- data.frame(
            bzx = bzx[filtered_index],
            bzy = bzy[filtered_index],
            bzx_se = bzx_se[filtered_index],
            bzy_se = bzy_se[filtered_index]
          )
          # 使用 bquote 函数创建动态标签
          x_label <- bquote(.(as.name(exp)) ~ (italic(b[zx])))
          y_label <- bquote(.(as.name(out)) ~ (italic(b[zy])))
          # 创建散点图
          p <- ggplot(data, aes(x = bzx, y = bzy)) +
            geom_point(shape = 20, size = 2, color = effect_col) +
            geom_segment(aes(x = bzx - bzx_se, xend = bzx + bzx_se, y = bzy, yend = bzy), size = 0.5, color = effect_col) +
            geom_segment(aes(x = bzx, xend = bzx, y = bzy - bzy_se, yend = bzy + bzy_se), size = 0.5, color = effect_col) +
            geom_abline(intercept = 0, slope = gsmr_results$bxy, linetype = "dashed", color = "dim grey") +
            labs(
              x = eval(x_label),
              y = eval(y_label)
            ) +
            annotate("text",
              x = min(data$bzx) * 1.2, y = max(data$bzy) * 1.2,
              label = paste(
                "OR:", ifelse(OR > 0.01, round(OR, 2), sprintf("%.2e", OR)),
                "\nP-value:", p_val_formatted
              ), hjust = 0, vjust = 1,
              size = 4
            )
          folder_name <- paste0(exp, "_to_", out)
          ggsave(filename = paste0(folder, "/GSMR_", folder_name, "_plot.pdf"), p, width = 8, height = 7)
        }
        gsmr_data <- data.frame(Exposure = exp, Outcome = out, b = b, se = se, OR = ORs, p = p_val_formatted)
        gsmr_data <- gsmr_data %>% dplyr::rename("OR(95%CI)" = OR, "P value" = p)
        gsmr_data_all <- rbind(gsmr_data_all, gsmr_data)
      }
    }
  }
  return(gsmr_data_all)
}
