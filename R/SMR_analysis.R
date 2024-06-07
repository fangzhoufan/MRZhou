#' Summary-data-based Mendelian Randomization(SMR) analysis
#'
#' @param smr_exe_dir the SMR path
#' @param g1000_dir the g1000/EUR path
#' @param qtls_dir the eqtl database,eg:eqtlgen,getxv8
#' @param input_dir the eqtl gene,eg:eqltgen:ENSG00000113161;GTEx:HMGCR
#' @param gwas_dir the vector of GWAS
#' @param outcome_prefix the file of results
#' @param multi_smr whether to analysis multi-smr
#' @param linux the performer choise
#'
#' @return the smr results
#' @export
#'
#' @examples
#' # SMR_analysis(qtls_dir='/public/home/fanfangzhou/R/GTEx8/SMR/eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense',gwas_dir=c('LDL.ma','RCS_UKB.ma','RCS_R10.ma'))
SMR_analysis <- function(smr_exe_dir = "/public/home/fanfangzhou/apps/SMR/smr_linux",
                         g1000_dir = "/public/home/fanfangzhou/R/1kgv3/g1000/EUR",
                         qtls_dir,
                         input_dir = "gene.list",
                         gwas_dir,
                         outcome_prefix = "targets_GWAS",
                         multi_smr = FALSE,
                         linux = TRUE) {
  input_dir <- file.path(getwd(), input_dir)
  output_dir <- file.path(getwd(), outcome_prefix, outcome_prefix)
  if (!dir.exists(output_dir)) {
    dir.create(outcome_prefix, recursive = TRUE)
  }

  # input_dir 这里，eqtlgen就要是如：ENSG00000113161，如果是GTEx就是gene symbol，如HMGCR
  besd_generate <- paste(
    shQuote(smr_exe_dir),
    "--beqtl-summary", shQuote(qtls_dir),
    "--query 1",
    "--extract-probe", shQuote(input_dir),
    "--out", shQuote(output_dir),
    "--make-besd"
  )
  # 使用system2执行chmod命令来添加执行权限
  if (linux) {
    # 写入到shell脚本文件
    script_file <- file.path(getwd(), outcome_prefix, "generate_besd.sh")
    writeLines(besd_generate, script_file)
    file_path <- paste0("'", script_file, "'")
    system2("chmod", args = c("+x", file_path), wait = TRUE)
    system2(script_file, wait = TRUE)
  } else {
    script_file <- file.path(getwd(), outcome_prefix, "generate_besd.bat")
    writeLines(besd_generate, script_file)
    file_path <- paste0("'", script_file, "'")
    system2(script_file, wait = TRUE)
  }


  cmd <- vector()
  # 循环处理每一个GWAS文件
  for (i in 1:length(gwas_dir)) {
    # 构建输出目录
    name <- gsub(".ma", "", gwas_dir[i])
    beqtl_dir <- file.path(getwd(), gwas_dir[i])
    output_dir <- file.path(getwd(), outcome_prefix, name)
    gwas_file <- file.path(getwd(), gwas_dir[i])
    # 构建完整的命令行
    cmd[i] <- paste(
      shQuote(smr_exe_dir),
      "--bfile", shQuote(g1000_dir),
      "--gwas-summary", shQuote(beqtl_dir),
      "--beqtl-summary", shQuote(output_dir),
      "--out", shQuote(output_dir)
    )
  }
  if (multi_smr) {
    cmd <- paste0(
      cmd,
      " --smr-multi",
      " --ld-multi-snp 0.1"
    )
  }

  if (linux) {
    # 写入到shell脚本文件
    script_file <- file.path(getwd(), outcome_prefix, "smr.sh")
    writeLines(cmd, script_file)
    file_path <- paste0("'", script_file, "'")
    system2("chmod", args = c("+x", file_path), wait = TRUE)
    system2(script_file, wait = TRUE)
  } else {
    script_file <- file.path(getwd(), outcome_prefix, "smr.bat")
    writeLines(besd_generate, script_file)
    file_path <- paste0("'", script_file, "'")
    system2(script_file, wait = TRUE)
  }


  # 检查命令是否成功执行
  if (result != 0) {
    stop("SMR analysis failed with error code: ", result, " for file: ", gwas_file)
  }

  message("SMR analysis completed successfully for all files.")
}
