#' PerformBinTestAIAnalysisForConditionNPointVect_knownCC
#'
#' Perform Binomial and QCC-corrected binomial tests with a given vector of point estimates, for given QCC.
#'
#' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
#' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
#' @param vectRepsCombsCC A vector of pairwise-computed correction constants for given replicates
#' @param ptVect A vector of point to compare with (should be compatible with the order and size of genes vector in table of allelic counts)
#' @param Q An optional parameter; quantile (for example 0.95, 0.8, etc)
#' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#' @param minDifference if specified, one additional column DAE is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to having non-overlapping CIs)
#' @return A table of gene names, AIs + CIs, classification into genes demonstrating signifficant difference (TRUE) from corresponding point estimate AI and those that don't (FALSE)
#' @export
#'
#' @importFrom stats "prop.test"
#'
PerformBinTestAIAnalysisForConditionNPointVect_knownCC <- function(inDF, vectReps,
                                                                   vectRepsCombsCC,
                                                                   ptVect,
                                                                   Q=0.95,
                                                                   thr=NA, thrUP=NA, thrType="each",
                                                                   minDifference=NA){

  DF <- ComputeAICIs(inDF, vectReps,
                     vectRepsCombsCC,
                     pt = ptVect,
                     Q=Q, BF=T,
                     thr=thr, thrUP=thrUP, thrType=thrType)

  DF$ptVect <- ptVect

  # Bonferroni correction:
  p_thr <- (1-Q)/nrow(na.omit(DF))

  # p-values -> tests:
  DF$BT <- (DF$BT_pval <= p_thr)
  DF$BT_CC <- (DF$BT_pval_CC <= p_thr)

  if (!is.na(minDifference))
  {
    DF$BT_CC_thrDiff <- (DF$BT_CC & (abs(DF$AI - ptVect) >= minDifference))
  }

  return(DF)
}
