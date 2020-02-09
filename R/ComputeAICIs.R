#' ComputeAICIs
#'
#' Calculate Binomial and QCC-corrected binomial CIs for a given vector of point estimates, for given QCC.
#'
#' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
#' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
#' @param vectRepsCombsCC A vector of pairwise-computed correction constants for given replicates
#' @param Q An optional parameter; quantile (for example 0.95, 0.8, etc)
#' @param BF Bonferroni correction, set False if Q is alredy corrected
#' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#' @return A table of gene names, AIs + CIs
#' @export
#'
#' @importFrom stats "prop.test"
#'
ComputeAICIs <- function(inDF, vectReps,
                         vectRepsCombsCC,
                         Q=0.95, BF=T,
                         thr=NA, thrUP=NA, thrType="each"){

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
  if (BF) {
    Qbf <- 1 - (1-Q)/nrow(na.omit(DF))
  } else {
    Qbf <- Q
  }

  # Bin test:
  tmpDFbt <- t(sapply(1:nrow(DF), function(i){
    if(is.na(matCOV[i]) | is.na(sumCOV[i]) | sumCOV[i]==0) { return(c(NA,NA)) }
    BT <- prop.test(matCOV[i], sumCOV[i], alternative="two.sided", conf.level = Qbf, correct=F)
    c(BT$conf.int[1], BT$conf.int[2])
  }))


  DF$BT_CIleft = tmpDFbt[, 1]
  DF$BT_CIright = tmpDFbt[, 2]

  DF$BT_CIleft_CC <- sapply(DF$AI - (DF$AI - tmpDFbt[, 1]) * CC,
                            function(lb){ max(0, lb) })
  DF$BT_CIright_CC <- sapply(DF$AI + (tmpDFbt[, 2] - DF$AI) * CC,
                             function(ub){ min(1, ub) })
  return(DF)
}
