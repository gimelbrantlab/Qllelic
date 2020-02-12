#' ComputeAICIs
#'
#' Calculate Binomial and QCC-corrected binomial CIs for a given vector of point estimates, for given QCC.
#'
#' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
#' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
#' @param vectRepsCombsCC A vector of pairwise-computed correction constants for given replicates
#' @param pt A value or a vector of values to compare with (if second, should be compatible with the order and size of genes vector in table of allelic counts)
#' @param Q An optional parameter; quantile (for example 0.95, 0.8, etc)
#' @param BF Bonferroni correction, set False if Q is alredy corrected
#' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#' @return A table of gene names, AI, coverage, test p-value, and Confidence Intervals
#' @export
#'
#' @importFrom stats "prop.test"
#'
ComputeAICIs <- function(inDF, vectReps,
                         vectRepsCombsCC,
                         pt = 0.5,
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

  # Point or points vector?

  edgeMover <- function(x){
    if(!is.na(x) & x == 0) {
      x = 0.00001
    } else if(!is.na(x) & x == 1) {
      x = 0.99999
    }
    return (x)
  }
  if (length(pt) == 1){
    pt_vect = rep(edgeMover(pt), times = nrow(DF))
  } else {
    pt_vect = sapply(pt, function(x){edgeMover(x)})
  }

  # Bin test:
  tmpDFbt <- t(sapply(1:nrow(DF), function(i){
    if(is.na(matCOV[i]) | is.na(sumCOV[i]) | sumCOV[i]==0) { return(c(NA,NA,NA)) }
    BT <- prop.test(matCOV[i], sumCOV[i], alternative="two.sided", p = pt_vect[i], conf.level = Qbf, correct=F)
    c(BT$p.value, BT$conf.int[1], BT$conf.int[2])
  }))

  DF$BT_CIleft  <- tmpDFbt[, 2]
  DF$BT_CIright <- tmpDFbt[, 3]

  # QCC-corrected Bin Test:
  k <- 1/CC**2

  tmpDFbtcc <- t(sapply(1:nrow(DF), function(i){
    if(is.na(matCOV[i]) | is.na(sumCOV[i]) | sumCOV[i]==0) { return(c(NA,NA,NA)) }
    BTcc <- prop.test(matCOV[i] * k, sumCOV[i] * k, alternative="two.sided", p = pt_vect[i], conf.level = Qbf, correct=F)
    c(BTcc$p.value, BTcc$conf.int[1], BTcc$conf.int[2])
  }))

  DF$BT_CIleft_CC  <- tmpDFbtcc[, 2]
  DF$BT_CIright_CC <- tmpDFbtcc[, 3]

  DF$BT_pval    <- tmpDFbt[, 1]
  DF$BT_pval_CC    <- tmpDFbtcc[, 1]

  return(DF)
}
