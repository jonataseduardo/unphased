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
#' population.
#'
#' @param prior_key string. Gives the which prior information should be used to
#' calculate the mutual information (genetic or uniform). The uniform prior
#' probaility for one population is equal to one divided by the number of
#' populations. The genetic prior probability is given by the number of samples
#' of that population divided by the total number of samples of all
#' populations. Default (prior_key = 'uniform').


popgen_mutual_information <- 
   function(freq_list, prior_key = 'uniform'){

    datum <- rbindlist(freq_list)
    datum[, N := .N, by = .(CHR, POS, SNP, VAR)]
    datum <- datum[N == length(freq_list)]

    if(prior_key == 'genetic'){
      datum[, qi := NCHROBS/sum(NCHROBS), 
            by = .(CHR, POS, SNP, VAR)]
    }
    if(prior_key == 'uniform'){
      datum[, qi := 1 / length(freq_list)]
    }

    out <- 
      datum[, .(pj = sum(AF * qi), 
                Ij = sum(AF * qi * log(AF))), 
             by = .(CHR, POS, SNP, VAR)
            ][, .(I1 = - sum(pj * log(pj)), I2 = sum(Ij)), 
             by = .(CHR, POS, SNP)
            ][, In := I1 + I2
            ][,`:=`(I1= NULL, I2 = NULL)]

    return(out[])
  }
