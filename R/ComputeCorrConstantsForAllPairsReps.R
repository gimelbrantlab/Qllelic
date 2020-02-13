#' ComputeCorrConstantsForAllPairsReps
#'
#' Computes QCC for all possible pairs of given replicates.
#'
#' @param inDF Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)
#' @param vectReps A vector of >= 2 replicate numbers for which the analysis should be applied
#' @param binNObs Optional (default=40), threshold on number of observations per bin
#' @param fitCovThr Optional (default=50), threshold on coverage for genes that will be included in Beta-Bin fitting
#' @param EPS Optional (default=1.05), base of exponent for the coverage binning, setting greater base (1.1 or 1.2 or 1.3) would result in fewer number of coverage bins in fitting process, thus will increase the computational speed, but may potentially reduce accuracy
#' @param thr Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis
#' @param thrUP Optional (default=NA), threshold for max gene coverage (default = NA)
#' @param thrType Optional (default = "each", also can be "average" for average coverage on replicates), threshold type

#' @return List of fitting outputs of ComputeCorrConstantFor2Reps() for each combinatorial pair of replicates (in order 1-2,1-3,..,1-N,2-3,..2-N,..,(N-1)-N): list with (1) fitted QCC ($fittedCC) and (2) a table with proportions of observed to expected quantiles per coverage bin ($QObsExpPropsTable).
#'
#' @export
#'
ComputeCorrConstantsForAllPairsReps <- function(inDF, vectReps, binNObs=40, fitCovThr=50,
                                                EPS=1.05, thr=NA, thrUP=NA, thrType="each"){

  options(stringsAsFactors = FALSE)

  repCombs <- combn(vectReps, 2)

  fitDATA <- lapply(1:ncol(repCombs), function(j){
    x = repCombs[, j]
    ComputeCorrConstantFor2Reps(inDF=inDF, reps=x, binNObs=binNObs, fitCovThr=fitCovThr,
                                EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType)
  })
  return(fitDATA)
}
