#' Moran I statitics
#'
#' \code{popphen_conditional_mutual_information()} Calculate genetic and population mutual
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
#' @param pheno_data data.table. Table with POP, mu, ssqr, NPHENO, columns
#' where POP has the population names, mu has the phenotypes means, ssqr has
#' the phenotypes standard deviations and NPHENO are the number of samples from
#' where the phenotypes where measured

moran_test <-
  function(data_in, col_name, wd_width = 200e3, alternative = "greater"){
    datum <- copy(data_in)
    setkey(datum, CHR, POS)

    datum[, aux := get(col_name)]
    datum[, xbar := mean(aux)]
    datum[, m2 := var((aux))]

    datum[, previous_POS := POS - wd_width / 2]
    datum[, next_POS := POS + wd_width / 2]
    datum[, idx := .I]
    datum[, z :=  (aux - xbar)]

    #stackoverflow.com/questions/32618803/
    #using-zoos-rollsum-within-data-table-on-timestamped-transactions

    resP <- datum[datum, .(idx = i.idx, seq = idx:i.idx),
                  by = .EACHI,
                  roll = -Inf, 
                  on = c(CHR = 'CHR', POS = 'previous_POS')]
     
    resN <- datum[datum, .(idx = i.idx, seq = idx:i.idx),
                  by = .EACHI,
                  roll = Inf, 
                  on = c(CHR = 'CHR', POS = 'next_POS')]
     
    res <- rbindlist(list(resN[idx != seq], resP))
    setkey(res, CHR, idx, seq)

    datum[, zw := datum[res$seq, 
                        sum(z)/ .N ^ 2, 
                        by = res$idx]$V1]

    datum[, n_neigh := datum[res$seq, 
                             .N, 
                             by = res$idx]$N]

    datum[, moran_I :=  z * zw / m2]

    if (alternative == "two.sided")
      datum[, p.value := - pnorm(moran_I, 
                                 0, 1, 
                                 low = FALSE, 
                                 log = TRUE) * log10(exp(1)) -log10(2)]
    else if (alternative == "greater")
      datum[, p.value := - pnorm(moran_I, 
                                 0, 1, 
                                 low = FALSE, 
                                 log = TRUE) * log10(exp(1))]
    else 
      datum[, p.value := - pnorm(-moran_I, 
                                 0, 1, 
                                 low = FALSE, 
                                 log = TRUE) * log10(exp(1))]

    datum[, c('xbar', 'm2', 
              'previous_POS', 'next_POS', 
              'z', 'zw', 'idx', col_name) := NULL]

    setnames(datum, 'aux', 'value')

    return(datum[,.(CHR, POS, SNP, value, n_neigh, moran_I, p.value)] )
  }
