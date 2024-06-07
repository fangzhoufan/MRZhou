#' SMR plot
#'
#' @param smr_exe_dir the SMR path
#' @param g1000_dir the g1000/EUR path
#' @param qtls_dir the eqtl database,eg:eqtlgen,getxv8
#' @param gwas_dir the vector of GWAS
#' @param queried_probe eg:ENSG00000113161
#' @param glist_dir eg:files of glist-hg19
#' @param outcome_prefix the file of results
#' @param gwas_label the plot title label
#' @param LocusPlot the choice of plot type
#' @param linux the performer choise
#' @param Rscript the Rcode of SMR plot
#'
#' @return a plot pdf
#' @export
#'
#' @examples
smr_plot <- function(smr_exe_dir = "/public/home/fanfangzhou/apps/SMR/smr_linux",
                     g1000_dir = "/public/home/fanfangzhou/R/1kgv3/g1000/EUR",
                     qtls_dir,
                     gwas_dir,
                     queried_probe,
                     glist_dir = "/public/home/fanfangzhou/R/GTEx8/SMR/glist-hg19",
                     outcome_prefix,
                     gwas_label, LocusPlot = F,
                     linux = TRUE,
                     Rscript = "/public/home/fanfangzhou/R/GTEx8/SMR/plot/plot/plot_SMR.r") {
  output_dir <- file.path(getwd(), outcome_prefix)
  if (!dir.exists(output_dir)) {
    dir.create(outcome_prefix, recursive = TRUE)
  }
  cmd <- vector()

  # 循环处理每一个GWAS文件
  for (i in 1:length(gwas_dir)) {
    # 构建输出目录
    name <- gsub(".ma", "", gwas_dir[i])
    gwas <- file.path(getwd(), gwas_dir[i])
    output_dir <- file.path(getwd(), outcome_prefix, name)
    # 构建完整的命令行
    cmd[i] <- paste(
      shQuote(smr_exe_dir),
      "--bfile", shQuote(g1000_dir),
      "--gwas-summary", shQuote(gwas),
      "--beqtl-summary", shQuote(qtls_dir),
      "--out", shQuote(output_dir),
      "--plot ",
      "--probe", shQuote(queried_probe),
      "--probe-wind 500 ",
      "--gene-list", shQuote(glist_dir)
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

  message("SMR analysis completed successfully for all files.")
  source(Rscript)
  for (i in 1:length(gwas_dir)){
    name<-gsub('.ma','',gwas_dir[i])
    filepath<-paste0(file.path(getwd(),outcome_prefix),'/plot/',name,'.',queried_probe,'.txt')
    SMRData = ReadSMRData(filepath)
    if(LocusPlot){
      pdf(paste0(file.path(getwd(),outcome_prefix),'/plot/',name,'.',queried_probe,'_SMRLocusPlot.pdf'),width =10,height =6)
      SMRLocusPlot(data=SMRData,
                   smr_thresh=8.4e-6, heidi_thresh=0.05, plotWindow=1000, max_anno_probe=16)
      dev.off()
    }else{
      pdf(paste0(file.path(getwd(),outcome_prefix),'/plot/',name,'.',queried_probe,'_SMREffectPlot.pdf'),width =7,height =6)
      SMREffectPlot(data=SMRData, trait_name=paste0('GWAS: ',gwas_label))
      dev.off()
    }
  }
}
