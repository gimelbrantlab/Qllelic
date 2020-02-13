#' PerformBinTestAIAnalysisForTwoConditions_knownCC
#'
#' Performs differential tests (with Bonferroni correction) for AI values for two conditions, for given QCC.
#'
#' @param inDF Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)
#' @param vect1CondReps A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
#' @param vect2CondReps A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
#' @param vect1CondRepsCombsCC A vector of pairwise-computed correction constants for first condition's tech reps (QCC=1 is no correction)
#' @param vect2CondRepsCombsCC A vector of pairwise-computed correction constants for second condition's tech reps (QCC=1 is no correction)
#' @param Q Optional (default=0.95), confidence level, quantile
#' @param thr Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis
#' @param thrUP Optional (default=NA), threshold for max gene coverage (default = NA)
#' @param thrType Optional (default = "each", also can be "average" for average coverage on replicates), threshold type
#' @param minDifference Optional (default=NA), if specified, one additional column is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to passing the test)
#'
#' @return A table of gene names, AIs + CIs for both conditions, p-values for both non-corrected (BT..) and QCC corrected (BT_CC..) differential tests, classification into genes demonstrating signifficant difference (TRUE) of AI estimates in two conditions, and those that don't (FALSE).
#'
#' @export
#'
#' @importFrom stats "prop.test"
#'
PerformBinTestAIAnalysisForTwoConditions_knownCC <- function(inDF, vect1CondReps, vect2CondReps,
                                                             vect1CondRepsCombsCC, vect2CondRepsCombsCC,
                                                             Q=0.95,
                                                             thr=NA, thrUP=NA, thrType="each", minDifference=NA){


  vectReps <- list(vect1CondReps, vect2CondReps)
  vectRepsCombsCC <- list(vect1CondRepsCombsCC,
                          vect2CondRepsCombsCC)
  meanRepsCombsCC <- c(mean(vect1CondRepsCombsCC), mean(vect2CondRepsCombsCC))

  DF_CI_divided <- lapply(1:2, function(i){
    ComputeAICIs(inDF, vectReps=vectReps[[i]], vectRepsCombsCC=vectRepsCombsCC[[i]],
                 Q=Q, BF=T, thr=thr, thrUP=thrUP, thrType=thrType)[, c("ID","sumCOV","matCOV","AI",
                                                                         "BT_CIleft","BT_CIright",
                                                                         "BT_CIleft_CC","BT_CIright_CC")]
  })

  DF <- merge(DF_CI_divided[[1]], DF_CI_divided[[2]], by = "ID")
  names(DF) <- c("ID",
                 paste0(names(DF_CI_divided[[1]])[-1], '_1'),
                 paste0(names(DF_CI_divided[[2]])[-1], '_2'))

  # Tests:
  diff_analysis_genes_num <- nrow(na.omit(DF))

  DF$BT_pval <- sapply(1:nrow(DF), function(i){
    if (!is.na(DF$matCOV_1[i]) & !is.na(DF$matCOV_2[i]) & DF$matCOV_1[i] > 0 & DF$matCOV_2[i] > 0){
      withCallingHandlers({
        prop.test(x = c(DF$matCOV_1[i], DF$matCOV_2[i]),
                  n = c(DF$sumCOV_1[i], DF$sumCOV_2[i]),
                  alternative="two.sided", correct=F)$p.value
      }, warning=function(w) {
        if (conditionMessage(w) == "Chi-squared approximation may be incorrect") {invokeRestart("muffleWarning")}
      })
    } else {
      NA
    }
  })

  k1 <- 1/meanRepsCombsCC[[1]]**2
  k2 <- 1/meanRepsCombsCC[[2]]**2

  DF$BT_CC_pval <- sapply(1:nrow(DF), function(i){
    if (!is.na(DF$matCOV_1[i]) & !is.na(DF$matCOV_2[i]) & DF$matCOV_1[i] > 0 & DF$matCOV_2[i] > 0){
      withCallingHandlers({
        prop.test(x = c(DF$matCOV_1[i] * k1, DF$matCOV_2[i] * k2),
                  n = c(DF$sumCOV_1[i] * k1, DF$sumCOV_2[i] * k2),
                  alternative="two.sided", correct=F)$p.value
      }, warning=function(w) {
        if (conditionMessage(w) == "Chi-squared approximation may be incorrect") {invokeRestart("muffleWarning")}
      })
    } else {
      NA
    }
  })

  DF$BT    <- (DF$BT_pval <= (1-Q)/diff_analysis_genes_num)
  DF$BT_CC <- (DF$BT_CC_pval <= (1-Q)/diff_analysis_genes_num)

  if (!is.na(minDifference))
  {
    DF$BT_CC_thrDiff <- (DF$BT_CC & (abs(DF$AI_1 - DF$AI_2) >= minDifference))
  }

  return(DF)
}
