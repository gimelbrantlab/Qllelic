#' PerformBinTestAIAnalysisForConditionNPoint_knownCC
#'
#' Perform Binomial and QCC-corrected binomial tests with a given point estimate, for given QCC.
#'
#' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
#' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
#' @param vectRepsCombsCC A vector of pairwise-computed correction constants for given replicates
#' @param pt A point to compare with
#' @param Q An optional parameter; quantile (for example 0.95, 0.8, etc)
#' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#' @param minDifference if specified, one additional column DAE is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to having non-overlapping CIs)
#' @return A table of gene names, AIs + CIs, classification into genes demonstrating signifficant difference from point estimate AI and those that don't
#' @export
#'
#' @importFrom stats "prop.test"
#'
PerformBinTestAIAnalysisForConditionNPoint_knownCC <- function(inDF, vectReps,
                                                               vectRepsCombsCC,
                                                               pt = 0.5,
                                                               Q=0.95,
                                                               thr=NA, thrUP=NA, thrType="each",
                                                               minDifference=NA){
  options(stringsAsFactors = FALSE)

  CC <- mean(vectRepsCombsCC)

  AI <- CountsToAI(inDF, reps=vectReps, meth="mergedToProportion", thr=thr, thrUP=thrUP, thrType=thrType)$AI
  tmpDF  <- ThresholdingCounts(inDF, reps=vectReps, thr=thr, thrUP=thrUP, thrType=thrType)
  sumCOV <- rowSums(tmpDF[, -1])
  if(ncol(tmpDF) == 3){
    matCOV <- tmpDF[, 2]
  } else {
    matCOV <- rowSums(tmpDF[, seq(2,ncol(tmpDF),2)])
  }

  DF <- data.frame(ID=inDF[,1], sumCOV=sumCOV, matCOV=matCOV, AI=AI)

  # Bonferroni correction:
  Qbf <- 1 - (1-Q)/nrow(na.omit(DF))

  if(pt == 0) {
    ptt = 0.00001
  } else if(pt == 1) {
    ptt = 0.99999
  } else {
    ptt = pt
  }
  # Bin test:
  tmpDFbt <- t(sapply(1:nrow(DF), function(i){
    if(is.na(matCOV[i]) | is.na(sumCOV[i]) | sumCOV[i]==0) { return(c(NA,NA,NA)) }
    BT <- prop.test(matCOV[i], sumCOV[i], alternative="two.sided", conf.level = Qbf, p=ptt, correct=F)
    c(BT$p.value, BT$conf.int[1], BT$conf.int[2])
  }))

  DF$BT_pval = tmpDFbt[, 1]
  DF$BT_CIleft = tmpDFbt[, 2]
  DF$BT_CIright = tmpDFbt[, 3]

  DF$BT_CIleft_CC <- sapply(DF$AI - (DF$AI - tmpDFbt[, 2]) * CC,
                            function(lb){ max(0, lb) })
  DF$BT_CIright_CC <- sapply(DF$AI + (tmpDFbt[, 3] - DF$AI) * CC,
                             function(ub){ min(1, ub) })

  # DF$BT <- (DF$BT_pval < (1-Q)/nrow(na.omit(DF)))
  DF$BT <- !(DF$BT_CIleft <= pt & DF$BT_CIright >= pt)
  # Find intersecting intervals > call them FALSE (non-rejected H_0)
  DF$BT_CC <- !(DF$BT_CIleft_CC <= pt & DF$BT_CIright_CC >= pt)

  if (!is.na(minDifference))
  {
    DF$BT_CC_thrDiff <- (DF$BT_CC & (abs(DF$AI - pt) >= minDifference))
  }

  return(DF)
}
