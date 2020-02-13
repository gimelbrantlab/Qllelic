#' PerformBinTestAIAnalysisForConditionNPoint_knownCC
#'
#' Performs Binomial and QCC-corrected binomial tests (with Bonferroni correction) with a given point estimate, for given QCC.
#'
#' @param inDF Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)
#' @param vectReps A vector (>=2) of replicate numbers for which the analysis should be applied
#' @param vectRepsCombsCC A vector of pairwise-computed correction constants for given replicates (QCC=1 is no correction)
#' @param pt Optional (default=0.5), a value to compare with
#' @param Q Optional (default=0.95), confidence level, quantile
#' @param thr Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis
#' @param thrUP Optional (default=NA), threshold for max gene coverage (default = NA)
#' @param thrType Optional (default = "each", also can be "average" for average coverage on replicates), threshold type
#' @param minDifference Optional (default=NA), if specified, one additional column is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to passing the test)
#'
#' @return A table of gene names, AIs + CIs, p-values for both non-corrected and (BT..) and QCC corrected (BT_CC..) tests, classification into genes demonstrating signifficant difference (TRUE) from point estimate AI and those that don't (FALSE).
#'
#' @export
#'
PerformBinTestAIAnalysisForConditionNPoint_knownCC <- function(inDF, vectReps,
                                                               vectRepsCombsCC,
                                                               pt = 0.5,
                                                               Q=0.95,
                                                               thr=NA, thrUP=NA, thrType="each",
                                                               minDifference=NA){

  DF <- ComputeAICIs(inDF, vectReps,
                     vectRepsCombsCC,
                     pt = pt,
                     Q=Q, BF=T,
                     thr=thr, thrUP=thrUP, thrType=thrType)

  # Bonferroni correction:
  p_thr <- (1-Q)/nrow(na.omit(DF))

  # p-values -> tests:
  DF$BT <- (DF$BT_pval <= p_thr)
  DF$BT_CC <- (DF$BT_pval_CC <= p_thr)

  if (!is.na(minDifference))
  {
    DF$BT_CC_thrDiff <- (DF$BT_CC & (abs(DF$AI - pt) >= minDifference))
  }

  return(DF)
}
