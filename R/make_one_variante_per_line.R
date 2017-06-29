#' Transform data from one SNP per line to one variant per line
#'
#' @import data.table
#'
#' @param freq_data data.table with one SNP per line and MAF information
#' @param inplace Boolean (FALSE by defaut)
#' 
#' @return data.table with one alelle per line with its frequency

make_one_allele_perline <- 
  function(freq_data, inplace = FALSE){

    if(inplace){
      f_data <- copy(freq_data)
    }else{
      f_data <- freq_data
    }

    f_data[, `:=`(AF1 = MAF, AF2 = 1 - MAF)]
    f_data[, MAF := NULL]

    data_long <- 
      melt(f_data, 
           measure = patterns("^A[12]", "^AF[12]"),
           value.name = c('VAR', 'AF'))
    data_long[, `:=`(variable = NULL)]
    return(data_long)
  }
