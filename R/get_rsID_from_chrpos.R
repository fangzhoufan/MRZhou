# The original version should do some manually working
# With more understandings about dbsnp (https://www.ncbi.nlm.nih.gov/snp/) and ensembl (http://rest.ensembl.org/), I re-wrote the function
# The function provides two options for database, one for dbsnp, another for ensembl. I prefer dbsnp and add additional function.
# Also, build is considered

#' How to get rsID from chromosome and position
#'
#' @param dat the dataframe of exposure or outcome
#' @param col_chr chromosome the default is "chr_col".
#' @param col_pos chromosome position
#' @param col_ref_allele Name of column with non effect allele. Must contain only the characters "A", "C", "T" or "G". The default is "refAllele".
#' @param col_alt_allele Name of column with effect allele. Must contain only the characters "A", "C", "T" or "G". The default is "altAllele"
#' @param col_snp Required name of column with SNP rs IDs. The default is "SNP".
#' @param build the version
#' @param pop investigated population
#' @param database dbsnp database
#'
#' @return a dataframe
#' @export
#'
#' @examples
get_rsID_from_chrpos <- function(dat,
                                 col_chr = "chr",
                                 col_pos = "start",
                                 col_ref_allele = "refAllele",
                                 col_alt_allele = "altAllele",
                                 col_snp = "SNP",
                                 build = "37",
                                 pop = "EUR",
                                 database = "dbsnp") {
  stopifnot(database %in% c("ensembl", "dbsnp"))
  library(tidyverse)
  if (database == "ensembl") {
    # Create and get a url
    server <- ifelse(build == "37",
      "http://grch37.rest.ensembl.org",
      "http://rest.ensembl.org"
    )

    query_term <- paste0(
      server, "/vep/human/region/", dat[[col_chr]], ":",
      dat[[col_pos]], "-", dat[[col_pos]], "/", dat[[col_alt_allele]], "?"
    )

    query_term_alt <- paste0(
      server, "/vep/human/region/", dat[[col_chr]], ":",
      dat[[col_pos]], "-", dat[[col_pos]], "/", dat[[col_ref_allele]], "?"
    )

    dat[[col_snp]] <- lapply(1:nrow(dat), function(i) {
      print(paste0("searching for No. ", i, " SNP"))
      query_res <- httr::GET(query_term[i], httr::content_type("application/json"))

      httr::warn_for_status(query_res)

      # Convert R objects from JSON
      query_res <- httr::content(query_res)
      res_df <- jsonlite::fromJSON(jsonlite::toJSON(query_res))
      snp <- res_df$colocated_variants[[1]][["id"]][[1]]
      if (is.null(snp)) {
        query_res <- httr::GET(query_term_alt[i], httr::content_type("application/json"))
        httr::warn_for_status(query_res)

        # Convert R objects from JSON
        query_res <- httr::content(query_res)
        res_df <- jsonlite::fromJSON(jsonlite::toJSON(query_res))
        snp <- res_df$colocated_variants[[1]][["id"]][[1]]
        if (is.null(snp)) {
          return(NA)
        } else {
          return(snp)
        }
      } else {
        return(snp)
      }

      # alleles <- unlist(str_split(res_df$allele_string[[1]],"/"))
      # ref_allele <- unlist(res_df$colocated_variants[[1]]$minor_allele)
      # alt_allele <- alleles[alleles != ref_allele]
      # alt_allele_freq <-res_df$colocated_variants[[1]][["frequencies"]][[alt_allele]][[str_to_lower("EUR")]][[1]]
    })

    dat[[col_snp]] <- unlist(dat[[col_snp]])
  }


  if (database == "dbsnp") {
    search_build <- ifelse(build == "37", "[POSITION_GRCH37]", "[Base Position]")

    query_term <- paste0(
      dat[[col_chr]], "[CHR] AND Homo[ORGN] AND ",
      dat[[col_pos]], search_build
    )

    SNP <- lapply(1:nrow(dat), function(i) {
      print(paste0("searching for No. ", i, " SNP"))
      snp <- unlist(rentrez::entrez_search(db = "snp", term = query_term[i])$ids)
      if (is.null(snp)) {
        return(NA)
      } else {
        return(paste0("rs", snp[length(snp)]))
      }
    })

    dat[[col_snp]] <- unlist(SNP)

    if (nrow(dat) != length(SNP[!is.na(SNP)])) {
      dat_with_snp <- dat[!is.na(dat[[col_snp]]), ]
      dat_without_snp <- dat[is.na(dat[[col_snp]]), ]
      query_term_alt <- paste0(
        dat_without_snp[[col_chr]], "[CHR] AND Homo[ORGN] AND ",
        dat_without_snp[[col_pos]] + 1, search_build
      )
      SNP <- lapply(1:nrow(dat_without_snp), function(i) {
        print(paste0("Researching for No. ", i, " SNP"))
        snp <- unlist(rentrez::entrez_search(db = "snp", term = query_term_alt[i])$ids)
        if (is.null(snp)) {
          return(NA)
        } else {
          return(paste0("rs", snp[length(snp)]))
        }
      })
      dat_without_snp[[col_snp]] <- unlist(SNP)
      dat <- rbind(dat_with_snp, dat_without_snp)
    }
  }

  dat
}
