
#' Calculate phenotipic and population mutual information conditional to the
#' genetic data.
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
#'
#' @param prior_key string. Gives the which prior information should be used to
#' calculate the mutual information (phenotipic, genetic or uniform). The
#' uniform prior probaility for one population is equal to one divided by the
#' number of The phenotipic/genetic prior probability is given by the number of
#' phenotipic/genetic samples of that population divided by the total number of
#' samples of all populations. Default (prior_key = 'uniform').

popphen_conditional_mutual_information <- 
   function(freq_list, pheno_data, prior_key = 'uniform'){

    datum <- 
      merge(rbindlist(freq_list), pheno_data,          
            by = c("POP"))

    NPOP <- datum[, .GRP, by = POP][, .N]

    if(prior_key == 'genetic'){
      datum[,qi := NCHROBS/sum(NCHROBS), 
            by = .(CHR, POS, SNP, VAR)]
    }
    if(prior_key == 'uniform'){
      datum[, qi := 1 / (NPOP)]
    }
    if(prior_key == 'phenotipic'){
      datum[,qi := NPHENO/sum(NPHENO), 
            by = .(CHR, POS, SNP, VAR)]
    }

    datum <- 
      datum[, N := .N, by = .(CHR, POS)
           ][N >= 2 * (NPOP)
           ][, N := NULL]

    setkey(datum, CHR, POS, VAR, POP)

    datum[, alpha := AF * qi]

    qdata <- 
      pheno_data[
        datum[, .(alpha = first(qi)), by = POP], on = .(POP)]
    
    info_integral <-
      function(datum, l_on){
        DT1 <- copy(datum)
        setkeyv(DT1, c(l_on, 'POP'))
        DT1[, idx := .GRP, by = POP]
        DT2 <- DT1[rep(1:.N, (NPOP))] 
        DT2[, idx := rep(1:(NPOP), each = (DT1[,.N]))]
        DT <- DT1[DT2, on = c(l_on, 'idx')]
        DT[, 
            {
              intg <- 
              function(x){
                (0.3989423 / ssqr) * 
                exp(- 0.5 * (x - mu) ^ 2 / (ssqr ^ 2)) * 
                log(0.3989423 * 
                    sum((i.alpha / i.ssqr) *
                         exp(- 0.5 * (x - i.mu) ^ 2 / (i.ssqr ^ 2)) 
                        )
                    )
              };
              integrate(Vectorize(intg), 
                        lower = mu - 6 *  ssqr, 
                        upper = mu + 6 *  ssqr)[1]
            }, 
            by = c(l_on, 'POP', 'mu', 'ssqr', 'alpha')][]
      }

    CI1 <- 
      info_integral(datum,
                    l_on = c('CHR', 'POS', 'SNP', 'VAR')
                    )[, .(I1 = - sum(alpha * value)), 
                       by = .(CHR, POS, SNP)]

    CI1[, I2 := qdata[, -0.5 * sum(alpha * ( 1 + log( 2 * pi * ssqr ^ 2 )))]]

    CI3 <- 
      datum[, .(pj = sum(alpha)), 
             by = .(CHR, POS, SNP, VAR)
            ][, .(I3 = sum(pj * log(pj))), 
             by = .(CHR, POS, SNP)]

    CI <- CI1[CI3, on = .(CHR, POS, SNP)]
    CI[, CIn := I1 + I2 + I3][]
    CI[, c('I1', 'I2', 'I3') := NULL]
    return(CI)
  }
