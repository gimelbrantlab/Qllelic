#' ComputeCorrConstantsForAllPairsReps
#'
#' Computes QCC for all possible pairs of given replicates.
#'
#' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
#' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
#' @param binNObs Threshold on number of observations per bin
#' @param fitCovThr Threshold on coverage for genes that will be included in Beta-Bin fitting
#' @param EPS An optional parameter to set a log window for coverage binning
#' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#' @return List of fitting outputs of ComputeCorrConstantFor2Reps() for each combinatorial pair of replicates (in order 1-2,1-3,..,1-N,2-3,..2-N,..,(N-1)-N): list with (1) fitted QCC ($fittedCC) and (2) a table with proportions of observed to expected quantiles per coverage bin ($QObsExpPropsTable).
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
