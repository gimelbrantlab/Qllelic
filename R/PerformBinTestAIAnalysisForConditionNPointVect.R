#' PerformBinTestAIAnalysisForConditionNPointVect
#'
#' Calculate QCC. Perform Binomial and QCC-corrected binomial tests with a given vector of point estimates.
#'
#' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
#' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
#' @param ptVect A vector of point to compare with (should be compatible with the order and size of genes vector in table of allelic counts)
#' @param binNObs Threshold on number of observations per bin
#' @param fitCovThr Threshold on coverage for genes that will be included in Beta-Bin fitting
#' @param Q An optional parameter; quantile (for example 0.95, 0.8, etc)
#' @param EPS An optional parameter to set a log window for coverage binning
#' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#' @param minDifference if specified, one additional column DAE is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to having non-overlapping CIs)
#' @return List of (1) fitted QCC for all combanatorial pairs of replicates ($CC), (2) ComputeCorrConstantsForAllPairsReps() output ($FitDATA), and (3) PerformBinTestAIAnalysisForConditionNPointVect_knownCC() output ($Output).
#' @export
#'
PerformBinTestAIAnalysisForConditionNPointVect <- function(inDF, vectReps,
                                                           ptVect,
                                                           binNObs=40, Q=0.95,
                                                           fitCovThr=50, EPS=1.3,
                                                           thr=NA, thrUP=NA, thrType="each",
                                                           minDifference=NA){

  fitDATA <- ComputeCorrConstantsForAllPairsReps(inDF, vectReps=vectReps, binNObs=binNObs, fitCovThr=fitCovThr,
                                                 EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType)

  vectRepsCombsCC <- sapply(fitDATA, function(fd){
    fd$fittedCC
  })

  RES <- PerformBinTestAIAnalysisForConditionNPointVect_knownCC(inDF, vectReps=vectReps,
                                                                vectRepsCombsCC=vectRepsCombsCC,
                                                                ptVect=ptVect, Q=Q,
                                                                thr=thr, thrUP=thrUP, thrType=thrType,
                                                                minDifference=minDifference)

  return(list(CC = vectRepsCombsCC,
              FitDATA = fitDATA,
              Output = RES))
}
