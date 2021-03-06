% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ComputeAICIs.R
\name{ComputeAICIs}
\alias{ComputeAICIs}
\title{ComputeAICIs}
\usage{
ComputeAICIs(
  inDF,
  vectReps,
  vectRepsCombsCC,
  pt = 0.5,
  Q = 0.95,
  BF = T,
  thr = NA,
  thrUP = NA,
  thrType = "each"
)
}
\arguments{
\item{inDF}{Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)}

\item{vectReps}{A vector of replicate numbers for which the analysis should be applied}

\item{vectRepsCombsCC}{A vector of pairwise-computed correction constants for given replicates (QCC=1 is no correction)}

\item{pt}{Optional (default=0.5), a value or a vector of values to compare with (if second, should be compatible with the order and size of genes vector in table of allelic counts)}

\item{Q}{Optional (default=0.95), confidence level, quantile}

\item{BF}{Optional (default=True), Bonferroni correction for multiple testing, set False ONLY IF Q is alredy corrected}

\item{thr}{Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis}

\item{thrUP}{Optional (default=NA), threshold for max gene coverage (default = NA)}

\item{thrType}{Optional (default = "each", also can be "average" for average coverage on replicates), threshold type}
}
\value{
A table with IDs, AI estimates, coverage, test p-value, and Confidence Intervals
}
\description{
Calculates Binomial and QCC-corrected binomial CIs for a given vector of AI estimates, and calculates test statistics for comparison with a point or vector of points, for given QCC.
}
