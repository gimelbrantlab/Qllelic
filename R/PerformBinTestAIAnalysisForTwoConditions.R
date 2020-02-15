#' PerformBinTestAIAnalysisForTwoConditions
#'
#' Calculates QCC. Performs differential tests (with Bonferroni correction) for AI values for two conditions.
#'
#' @param inDF Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)
#' @param vect1CondReps A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
#' @param vect2CondReps A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
#' @param binNObs Optional (default=40), threshold on number of observations per bin
#' @param fitCovThr Optional (default=50), threshold on coverage for genes that will be included in Beta-Bin fitting
#' @param Q Optional (default=0.95), confidence level, quantile
#' @param EPS Optional (default=1.05), base of exponent for the coverage binning, setting greater base (1.1 or 1.2 or 1.3) would result in fewer number of coverage bins in fitting process, thus will increase the computational speed, but may potentially reduce accuracy
#' @param thr Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis
#' @param thrUP Optional (default=NA), threshold for max gene coverage (default = NA)
#' @param thrType Optional (default = "each", also can be "average" for average coverage on replicates), threshold type
#' @param minDifference Optional (default=NA), if specified, one additional column is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to passing the test)
#'
#' @return List of (1) fitted QCC for all combanatorial pairs of replicates for both conditions ($CC), (2) ComputeCorrConstantsForAllPairsReps() output for both conditions ($FitDATA), and (3) PerformBinTestAIAnalysisForTwoConditions_knownCC() output ($Output).
#'
#' @export
#'
PerformBinTestAIAnalysisForTwoConditions <- function(inDF, vect1CondReps, vect2CondReps, binNObs=40, fitCovThr=50, Q=0.95, EPS=1.05,
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
  #print(paste(vect1CondRepsCombsCC, vect2CondRepsCombsCC))

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
