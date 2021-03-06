% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maf_filter.R
\name{maf_filter}
\alias{maf_filter}
\title{Filter from data set SNPs with MAF bellow a threshold}
\usage{
maf_filter(freq_data, maf_th = 0.05, inplace = FALSE)
}
\arguments{
\item{freq_data}{data.table. Data with one allele per line.}

\item{maf_th}{float. MAF threshold. Default maf_th = 0.05}
}
\value{
data.table
}
\description{
\code{maf_filter()} Filter SNPs with MAF bellow a threshold
}
