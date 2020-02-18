#' ComputeAICIs
#'
#' Calculates Binomial and QCC-corrected binomial CIs for a given vector of AI estimates, and calculates test statistics for comparison with a point or vector of points, for given QCC.
#'
#' @param inDF Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)
#' @param vectReps A vector of replicate numbers for which the analysis should be applied
#' @param vectRepsCombsCC A vector of pairwise-computed correction constants for given replicates (QCC=1 is no correction)
#' @param pt Optional (default=0.5), a value or a vector of values to compare with (if second, should be compatible with the order and size of genes vector in table of allelic counts)
#' @param Q Optional (default=0.95), confidence level, quantile
#' @param BF Optional (default=True), Bonferroni correction for multiple testing, set False ONLY IF Q is alredy corrected
#' @param thr Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis
#' @param thrUP Optional (default=NA), threshold for max gene coverage (default = NA)
#' @param thrType Optional (default = "each", also can be "average" for average coverage on replicates), threshold type
#'
#' @return A table with IDs, AI estimates, coverage, test p-value, and Confidence Intervals
#'
#' @export
#'
#' @importFrom stats "prop.test"
#'
ComputeAICIs <- function(inDF, vectReps, vectRepsCombsCC, pt = 0.5, Q=0.95, BF=T, thr=NA, thrUP=NA, thrType="each"){

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
    BT <- withCallingHandlers({
        prop.test(matCOV[i], sumCOV[i], alternative="two.sided", p = pt_vect[i], conf.level = Qbf, correct=F)
      }, warning=function(w) {
        if (conditionMessage(w) == "Chi-squared approximation may be incorrect") {invokeRestart("muffleWarning")}
      })
    c(BT$p.value, BT$conf.int[1], BT$conf.int[2])
  }))

  DF$BT_CIleft  <- tmpDFbt[, 2]
  DF$BT_CIright <- tmpDFbt[, 3]

  # QCC-corrected Bin Test:
  k <- 1/CC**2

  tmpDFbtcc <- t(sapply(1:nrow(DF), function(i){
    if(is.na(matCOV[i]) | is.na(sumCOV[i]) | sumCOV[i]==0) { return(c(NA,NA,NA)) }
    BTcc <- withCallingHandlers({
      prop.test(matCOV[i] * k, sumCOV[i] * k, alternative="two.sided", p = pt_vect[i], conf.level = Qbf, correct=F)
    }, warning=function(w) {
      if (conditionMessage(w) == "Chi-squared approximation may be incorrect") {invokeRestart("muffleWarning")}
    })
    c(BTcc$p.value, BTcc$conf.int[1], BTcc$conf.int[2])
  }))

  DF$BT_CIleft_CC  <- tmpDFbtcc[, 2]
  DF$BT_CIright_CC <- tmpDFbtcc[, 3]

  DF$BT_pval    <- tmpDFbt[, 1]
  DF$BT_pval_CC    <- tmpDFbtcc[, 1]

  return(DF)
}
