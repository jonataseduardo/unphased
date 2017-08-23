#' Create data.table with CHR, SNP, POS, AF, POP columns
#'
#'\code{create_data_from_plink_files()} reads .bim and .frq plink files
#' and merge then into one data.table
#' 
#' @import data.table
#' @export
#' 
#' @param files_path string with full path of .bim and .frq files 
#' @param pop_name string with the name of the population. DEFAULT is the
#' prefix of .bim and .frq files.  
#' @param one_allele_perline boolean. Transform data from one SNP per line to
#' one allele per line, Default TRUE.
#' 
#' @return data.table with the information of .bim and .frq files 

create_data_from_plink_files <- 
  function(files_path, 
           pop_name = tail(unlist(strsplit(files_path, '/')), n = 1),
           one_allele_perline = TRUE){

    freq_data  <- 
      merge(
        fread(paste0(files_path, ".bim"),
              col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2")),
        fread(paste0(files_path, ".frq")), 
        by = c("CHR", "SNP", "A1", "A2"))

    freq_data[, `:=`(POP, pop_name)]
    if(one_allele_perline)
      freq_data <- make_one_allele_perline(freq_data)

    setkey(freq_data, CHR, POS)

    return(freq_data[])
  }
