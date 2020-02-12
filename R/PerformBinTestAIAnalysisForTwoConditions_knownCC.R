#' PerformBinTestAIAnalysisForTwoConditions_knownCC
#'
#' Perform differential tests for AI values for two conditions, for given QCC.
#'
#' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
#' @param vect1CondReps A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
#' @param vect2CondReps A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
#' @param vect1CondRepsCombsCC A vector of pairwise-computed correction constants for first condition's tech reps
#' @param vect2CondRepsCombsCC A vector of pairwise-computed correction constants for second condition's tech reps
#' @param Q An optional parameter; quantile (for example 0.95, 0.8, etc)
#' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#' @param minDifference if specified, one additional column DAE is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to having non-overlapping CIs)
#' @return A table of gene names, AIs + CIs for both conditions, classification into genes demonstrating signifficant difference of AI estimates in two conditions, and those that don't
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
      prop.test(x = c(DF$matCOV_1[i], DF$matCOV_2[i]),
                n = c(DF$sumCOV_1[i], DF$sumCOV_2[i]),
                alternative="two.sided", correct=F)$p.value
    } else {
      NA
    }
  })

  k1 <- 1/meanRepsCombsCC[[1]]**2
  k2 <- 1/meanRepsCombsCC[[2]]**2

  DF$BT_CC_pval <- sapply(1:nrow(DF), function(i){
    if (!is.na(DF$matCOV_1[i]) & !is.na(DF$matCOV_2[i]) & DF$matCOV_1[i] > 0 & DF$matCOV_2[i] > 0){
      prop.test(x = c(DF$matCOV_1[i] * k1, DF$matCOV_2[i] * k2),
                n = c(DF$sumCOV_1[i] * k1, DF$sumCOV_2[i] * k2),
                alternative="two.sided", correct=F)$p.value
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
