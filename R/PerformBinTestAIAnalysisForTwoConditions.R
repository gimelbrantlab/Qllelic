#' PerformBinTestAIAnalysisForTwoConditions
#'
#' Calculate QCC. Perform differential tests for AI values for two conditions.
#'
#' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
#' @param vect1CondReps A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
#' @param vect2CondReps A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
#' @param binNObs Threshold on number of observations per bin
#' @param fitCovThr Threshold on coverage for genes that will be included in Beta-Bin fitting
#' @param Q An optional parameter; quantile (for example 0.95, 0.8, etc)
#' @param EPS An optional parameter to set a log window for coverage binning
#' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#' @param minDifference if specified, one additional column DAE is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to having non-overlapping CIs)
#' @return List of (1) fitted QCC for all combanatorial pairs of replicates for both conditions ($CC), (2) ComputeCorrConstantsForAllPairsReps() output for both conditions ($FitDATA), and (3) PerformBinTestAIAnalysisForTwoConditions_knownCC() output ($Output).
#' @export
#'
PerformBinTestAIAnalysisForTwoConditions <- function(inDF, vect1CondReps, vect2CondReps, binNObs=40, fitCovThr=50, Q=0.95, EPS=1.3,
                                                     thr=NA, thrUP=NA, thrType="each", minDifference=NA){

  fitDATA1Cond <- ComputeCorrConstantsForAllPairsReps(inDF, vectReps=vect1CondReps, binNObs=binNObs, fitCovThr=fitCovThr,
                                                      EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType)
  fitDATA2Cond <- ComputeCorrConstantsForAllPairsReps(inDF, vectReps=vect2CondReps, binNObs=binNObs, fitCovThr=fitCovThr,
                                                      EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType)

  vect1CondRepsCombsCC <- sapply(fitDATA1Cond, function(fd){
    fd$fittedCC
  })
  vect2CondRepsCombsCC <- sapply(fitDATA2Cond, function(fd){
    fd$fittedCC
  })
  print(paste(vect1CondRepsCombsCC, vect2CondRepsCombsCC))

  RES <- PerformBinTestAIAnalysisForTwoConditions_knownCC(inDF, vect1CondReps=vect1CondReps, vect2CondReps=vect2CondReps,
                                                          vect1CondRepsCombsCC=vect1CondRepsCombsCC,
                                                          vect2CondRepsCombsCC=vect2CondRepsCombsCC,
                                                          Q=Q,
                                                          thr=thr, thrUP=thrUP, thrType=thrType,
                                                          minDifference=minDifference)

  return(list(CC = list(vect1CondRepsCombsCC, vect2CondRepsCombsCC),
              FitDATA = list(fitDATA1Cond, fitDATA2Cond),
              Output = RES))
}
