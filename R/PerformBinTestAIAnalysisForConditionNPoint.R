#' PerformBinTestAIAnalysisForConditionNPoint
#'
#' Calculates QCC. Performs Binomial and QCC-corrected binomial tests (with Bonferroni correction) with a given point estimate.
#'
#' @param inDF Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)
#' @param vectReps A vector (>=2) of replicate numbers for which the analysis should be applied
#' @param pt Optional (default=0.5), a value to compare with
#' @param binNObs Optional (default=40), threshold on number of observations per bin
#' @param fitCovThr Optional (default=50), threshold on coverage for genes that will be included in Beta-Bin fitting
#' @param Q Optional (default=0.95), confidence level, quantile
#' @param EPS Optional (default=1.05), base of exponent for the coverage binning, setting greater base (1.1 or 1.2 or 1.3) would result in fewer number of coverage bins in fitting process, thus will increase the computational speed, but may potentially reduce accuracy
#' @param thr Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis
#' @param thrUP Optional (default=NA), threshold for max gene coverage (default = NA)
#' @param thrType Optional (default = "each", also can be "average" for average coverage on replicates), threshold type
#' @param minDifference Optional (default=NA), if specified, one additional column is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to passing the test)
#'
#' @return List of (1) fitted QCC for all combanatorial pairs of replicates ($CC), (2) ComputeCorrConstantsForAllPairsReps() output ($FitDATA), and (3) PerformBinTestAIAnalysisForConditionNPoint_knownCC() output ($Output).
#'
#' @export
#'
PerformBinTestAIAnalysisForConditionNPoint <- function(inDF, vectReps,
                                                       pt = 0.5,
                                                       binNObs=40,
                                                       Q=0.95,
                                                       fitCovThr=50, EPS=1.05,
                                                       thr=NA, thrUP=NA, thrType="each",
                                                       minDifference=NA){

  fitDATA <- ComputeCorrConstantsForAllPairsReps(inDF, vectReps=vectReps, binNObs=binNObs, fitCovThr=fitCovThr,
                                                 EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType)

  vectRepsCombsCC <- sapply(fitDATA, function(fd){
    fd$fittedCC
  })

  RES <- PerformBinTestAIAnalysisForConditionNPoint_knownCC(inDF, vectReps=vectReps,
                                                            vectRepsCombsCC=vectRepsCombsCC,
                                                            pt=pt, Q=Q,
                                                            thr=thr, thrUP=thrUP, thrType=thrType,
                                                            minDifference=minDifference)

  return(list(CC = vectRepsCombsCC,
              FitDATA = fitDATA,
              Output = RES))
}
