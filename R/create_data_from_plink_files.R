#' Create data.table with CHR, SNP, POS, AF, POP columns
#'
#'\code{create_data_from_plink_files()} reads .bim and .frq plink files
#' and merge then into one data.table
#' 
#' @param files_path string with full path of .bim and .frq files 
#' @param pop_name string with the name of the population. DEFAULT is the
#' prefix of .bim and .frq files  
#' 
#' @return data.table with the information of .bim and .frq files 
#' 
#' @import data.table
#' @export

create_data_from_plink_files <- 
  function(files_path, 
           pop_name = tail(unlist(strsplit(files_path, '/')), n = 1)){

    freq_data  <- 
      merge(
        fread(paste0(files_path, ".bim"),
              col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2")),
        fread(paste0(files_path, ".frq")), 
        by = c("CHR", "SNP", "A1", "A2"))

    freq_data[, `:=`(POP, pop_name)]
    return(freq_data)
  }
