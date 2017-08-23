#' Calculate genetic and population mutual information.
#'
#' \code{popgen_mutual_information()} Calculate genetic and population mutual
#' information given the allele frequecies and prior probability of an idividual
#' to belong to one population 
#'
#'
#' @import data.table
#' @export 
#' 
#' @param freq_list list of data.table. Each element of the list is a
#' data.table object with the allele frequency of the variants of one
#' population. The populations in the list should be in the order (focal,
#' ancestral, distal)

info_length <-
  function(data_list){

      info_123 <-
        popgen_mutual_information(data_list)
      
      info_23 <- 
        popgen_mutual_information(data_list[2:3])

      info_Pin <- 
        info_123[info_23, 
                 PIn := In  - i.In,
                 on = .(CHR, POS, SNP)][
                 !is.na(PIn)][, In := NULL][]
      return(info_PIn)
  }
