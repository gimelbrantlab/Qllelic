#
# PERFORM DIFF AI ANALYSIS ON 2 CONDITIONS or CONDITION AND POINT
# _______________________________________________________________________________________


# _______________________________________________________________________________________

options(stringsAsFactors = FALSE)
# _______________________________________________________________________________________

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: BETA-BINOMIAL FITTING
# ---------------------------------------------------------------------------------------

# wi -- weight of betabin(alphai); wi>0; w1+w2=1
# alphai -- a=b coefficient in beta
# c -- coverage
# N -- number of observations (n)
# K -- (=2) number of classes (k)

LogFraction <- function(a_up, a_down, j){
  #' Input: 3 numbers (double, double, int)
  #'
  #' @param a_up A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param a_down A vector (2) of replicate numbers that should be considered
  #' @param j
  #' @return log((j + a_up) / (j + a_down))
  #' @examples
  #'
  return(log1p((a_up -1)/(j+1)) - log1p((a_down -1)/(j+1)))
}

MixBetaBinomialFitStep <- function(initials_old, coverage, observations){
  #' Input: 3 numbers (double, double, int)
  #'
  #' @param initials_old Initials for EM step: initials = c(w1, alpha1, alpha2), weight of first component and alphas for both beta-binomial distributions in a mixture
  #' @param coverage A number, that reflects higher boundary of the coerage bin
  #' @param observations A vector of maternal counts in the bin
  #' @return Re-fitted initials for next EM step.
  #' @examples
  #'
  w <- c(initials_old[1], 1-initials_old[1])
  alpha <- initials_old[2:3]

  # EM step 1:
  gamma_n_k <- sapply(1:2, function(k){
    a_k <- alpha[k]
    a_notk <- alpha[2-(k+1)%%2]
    w_k <- w[k]
    w_notk <- w[2-(k+1)%%2]

    gamma_k <- sapply(observations, function(xn){
      if(xn == 0){
        logFractionsProd <- sum(sapply(0:(coverage-xn-1), function(j){LogFraction(a_notk, a_k, j)})) +
          sum(sapply(0:(coverage-1), function(j){LogFraction(2*a_k, 2*a_notk, j)}))
      } else if(xn == coverage) {
        logFractionsProd <- sum(sapply(0:(xn-1), function(j){LogFraction(a_notk, a_k, j)})) +
          sum(sapply(0:(coverage-1), function(j){LogFraction(2*a_k, 2*a_notk, j)}))
      } else {
        logFractionsProd <- sum(sapply(0:(xn-1), function(j){LogFraction(a_notk, a_k, j)})) +
          sum(sapply(0:(coverage-xn-1), function(j){LogFraction(a_notk, a_k, j)})) +
          sum(sapply(0:(coverage-1), function(j){LogFraction(2*a_k, 2*a_notk, j)}))
      }
      result <- 1 / (1 + w_notk/w_k * exp(logFractionsProd))
      ########################### 1 + eps != 1 (!), ну и чёрт с ним.
      return(result)
    })
    return(gamma_k)
  })

  # EM step 2:
  sum_gamma_k <- colSums(gamma_n_k)

  w_new <- sum_gamma_k / length(observations)
  #mean_new <- sapply(1:2, function(k){
  #  sum(gamma_n_k[, k] * observations) / sum_gamma_k[k]
  #})
  mean_new <- c(0.5*coverage, 0.5*coverage)
  var_new <- sapply(1:2, function(k){
    sum(gamma_n_k[, k] * (observations - mean_new[k])**2) / sum_gamma_k[k]
  })
  alpha_new <- sapply(1:2, function(k){
    (coverage**2 - 4 * var_new[k]) / (8 * var_new[k] - 2 * coverage)
  })

  initials_new <- c(w_new[1], alpha_new)
  return(as.numeric(initials_new))
}

MixBetaBinomialFit <- function(initials, coverage, observations){
  #' Input: 3 numbers (double, double, int)
  #'
  #' @param initials Initials for EM algm: initials = c(w1, alpha1, alpha2), weight of first component and alphas for both beta-binomial distributions in a mixture
  #' @param coverage A number, that reflects higher boundary of the coerage bin
  #' @param observations A vector of maternal counts in the bin
  #' @return Weight of first component and alphas for both beta-binomial distributions in a mixture, to which algm coincided, plus number of steps.
  #' @examples
  #'

  initials_old <- initials
  initials_new <- MixBetaBinomialFitStep(initials_old, coverage, observations)
  n_steps <- 1
  while(all(abs(initials_new - initials_old) > 0.001) &
        initials_new[3] < 2){
    initials_old <- initials_new
    initials_new <- MixBetaBinomialFitStep(initials_old, coverage, observations)
    n_steps = n_steps + 1
  }
  return(c(initials_new, n_steps))
}


# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: COMPUTE CORR CONSTANT
# ---------------------------------------------------------------------------------------

ComputeCorrConstantFor2Reps <- function(inDF, reps, binNObs=40,
                                        EPS=1.3, thr=NA, thrUP=NA, thrType="each"){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for each condition
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param reps A vector of 2 of replicate numbers that should be considered
  #' @param binNObs Threshold on number of observations per bin
  #' @param EPS An optional parameter to set a log window for coverage binning
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return A table of gene names, AIs + CIs for each condition, classification into genes demonstrating differential AI and those that don't
  #' @examples
  #'

  ##-------------------------------------------------------------------------------------------------------------------------------------
  ## 1. Compute beta-binomial parameters for merged AI distribution per each coverage bin:
  ##-------------------------------------------------------------------------------------------------------------------------------------
  ##     1.1. Table with counts and ai:
  ##

  meancov =  MeanCoverage(inDF, reps=reps, thr=thr, thrUP=thrUP, thrType=thrType)$meanCOV
  df_unit_info = data.frame(ID = inDF[, 1],
                            AI = CountsToAI(inDF, reps=reps, thr=thr, thrUP=thrUP, thrType=thrType)$AI,
                            AI_merged = CountsToAI(inDF, reps=reps, meth="mergedToProportion", thr=thr, thrUP=thrUP, thrType=thrType)$AI,
                            AI_1 = CountsToAI(inDF, reps=reps[1], thr=thr, thrUP=thrUP, thrType=thrType)$AI,
                            AI_2 = CountsToAI(inDF, reps=reps[2], thr=thr, thrUP=thrUP, thrType=thrType)$AI,
                            mCOV = meancov,
                            COV = meancov*2,
                            binCOV = ceiling(EPS**ceiling(log(meancov, base=EPS))))
  df_covbinsnum = data.frame(binCOV = sort(unique(df_unit_info$binCOV)),
                             binNUM = sapply(sort(unique(df_unit_info$binCOV)), function(x){
                               sum(!is.na(df_unit_info$binCOV) & df_unit_info$binCOV == x)
                             })
  )
  df_unit_info = na.omit(df_unit_info)
  df_unit_info$binNUM = sapply(df_unit_info$binCOV, function(x){
    df_covbinsnum[df_covbinsnum$binCOV==x, ]$binNUM
  })

  ##     1.2. Calculate parameters:
  ##

  initials = c(0.5, c(10, 1/50))
  if (is.na(thr)) {
    thr = 0
  }

  covbinsGthr = df_covbinsnum$binCOV[df_covbinsnum$binNUM > binNObs &
                                       df_covbinsnum$binCOV >= max(50, thr)]

  print(paste(length(covbinsGthr), "COVERAGE BINS"))

  df_betabin_params = do.call(rbind, lapply(1:length(covbinsGthr), function(i){
    coverage = covbinsGthr[i]
    df = df_unit_info[df_unit_info$binCOV == coverage, ]
    observations = round(df$AI_merged * coverage*2)

    fitres = MixBetaBinomialFit(initials, coverage*2, observations)
    dfres = data.frame(No = i,
                       coverage = coverage,
                       x_weight_predicted = fitres[1],
                       x_alpha_predicted = fitres[2],
                       y_alpha_predicted = fitres[3],
                       algm_steps = fitres[4])
    print(paste(i, "|", "meanCOV:", covbinsGthr[i], ",", "#STEPS:", dfres$algm_steps, sep="    "))
    return(dfres)
  }))

  ##-------------------------------------------------------------------------------------------------------------------------------------
  ## 2. Simulate AI distribution based on beta patameters for each bin:
  ##-------------------------------------------------------------------------------------------------------------------------------------
  ##     2.1. Generating beta-distributed probabilities for each coverage bin:
  ##

  lst_ai_betafit = lapply(1:nrow(df_betabin_params), function(i){
    A = rbeta(5000*df_betabin_params$x_weight_predicted[i],
              df_betabin_params$x_alpha_predicted[i],
              df_betabin_params$x_alpha_predicted[i])
    B = rbeta(5000*(1 - df_betabin_params$x_weight_predicted[i]),
              df_betabin_params$y_alpha_predicted[i],
              df_betabin_params$y_alpha_predicted[i])
    c(A,B)
  })

  ##     2.2. Generating pairs of betabin-distributed AIs for each coverage bin and prob-s:
  ##

  lst_2ai_betabinfit = lapply(1:length(covbinsGthr), function(i){
    sapply(lst_ai_betafit[[i]], function(p){
      rbinom(2, covbinsGthr[i], p) / covbinsGthr[i]
    })
  })

  ##     2.3. Compute differences for fitted AI for each coverage bin:
  ##

  lst_dai_betabinfit = lapply(lst_2ai_betabinfit, function(df){
    abs(df[1, ] - df[2, ])
  })

  ##      2.4 Compute real AI differences for each coverage bin:
  ##

  lst_dai_observed = lapply(covbinsGthr, function(covbin){
    df = df_unit_info[df_unit_info$binCOV == covbin, ]
    abs(df$AI_1 - df$AI_2)
  })


  ##-------------------------------------------------------------------------------------------------------------------------------------
  ## 3. Compute Q% quantiles for both sets and a ratio (CC estimates per bin and Q):
  ##-------------------------------------------------------------------------------------------------------------------------------------

  QQ = c(0.2,0.35,0.5,0.65,0.8,0.9,0.95)

  df_observed_expected_quantile_proportions = do.call(rbind, lapply(1:length(covbinsGthr), function(c){
    do.call(rbind, lapply(1:length(QQ), function(q){
      df = data.frame(Q = QQ[q],
                      binCOV = covbinsGthr[c],
                      QV_obs = quantile(lst_dai_observed[[c]], QQ[q], na.rm = T),
                      QV_fit = quantile(lst_dai_betabinfit[[c]], QQ[q], na.rm = T))
      df$CC = df$QV_obs / df$QV_fit
      df
    }))
  }))

  #print(df_observed_expected_quantile_proportions)

  ##-------------------------------------------------------------------------------------------------------------------------------------
  ## 4. Fit the ratio observed/predicted:
  ##-------------------------------------------------------------------------------------------------------------------------------------

  df_observed_expected_quantile_proportions$CC[!is.finite(df_observed_expected_quantile_proportions$CC)] = NA

  fittedCC = lm(data = df_observed_expected_quantile_proportions,
                CC ~ 1,
                na.action=na.exclude)$coefficients[1]

  print(fittedCC)
  ##-------------------------------------------------------------------------------------------------------------------------------------
  return(list(
    infoTable = df_unit_info,
    betaBinTable = df_betabin_params,
    dAIBetaBinFit = lst_dai_betabinfit,
    dAIObserved = lst_dai_observed,
    QObsExpPropsTable = df_observed_expected_quantile_proportions,
    fittedCC = fittedCC
  ))

}


ComputeCorrConstantsForAllPairsReps <- function(inDF, vectReps, binNObs=40,
                                                EPS=1.3, thr=NA, thrUP=NA, thrType="each"){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for each condition
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param binNObs Threshold on number of observations per bin
  #' @param EPS An optional parameter to set a log window for coverage binning
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return A table of gene names, AIs + CIs for each condition, classification into genes demonstrating differential AI and those that don't
  #' @examples
  #'

  repCombs <- combn(vectReps, 2)

  fitDATA <- lapply(1:ncol(repCombs), function(j){
    x = repCombs[, j]
    ComputeCorrConstantFor2Reps(inDF=inDF, reps=x, binNObs=binNObs,
                                EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType)
  })
  return(fitDATA)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PERFORM DIFF CI(AI) ANALYSIS -- CONDITION & POINT
# ---------------------------------------------------------------------------------------

PerformBinTestAIAnalysisForConditionNPoint_knownCC <- function(inDF, vectReps, vectRepsCombsCC,
                                                               pt = 0.5, binNObs=40, Q=0.95,
                                                               thr=NA, thrUP=NA, thrType="each", minDifference=NA){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for condition + point estimate to compare
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param vectRepsCombsCC A vector of pairwise-computed correction constants for given replicates
  #' @param pt A point to compare with
  #' @param binNObs Threshold on number of observations per bin
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @param minDifference if specified, one additional column DAE is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to having non-overlapping CIs)
  #' @return A table of gene names, AIs + CIs, classification into genes demonstrating differential from point estimate AI and those that don't
  #' @examples
  #'

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

  # Bin test:
  tmpDFbt <- t(sapply(1:nrow(DF), function(i){
    if(is.na(matCOV[i]) | is.na(sumCOV[i]) | sumCOV[i]==0) { return(c(NA,NA,NA)) }
    BT <- binom.test(matCOV[i], sumCOV[i], alternative="two.sided", conf.level = Qbf, p=pt)
    c(BT$p.value, BT$conf.int[1], BT$conf.int[2])
  }))

  # DF$BT_pval = tmpDFbt[, 1]
  # DF$BT_CIleft = DF$AI - (DF$AI - tmpDFbt[, 2]) / sqrt(length(vectReps))
  # DF$BT_CIright = DF$AI + (tmpDFbt[, 3] - DF$AI) / sqrt(length(vectReps))
  #
  # DF$BT_CIleft_CC <- sapply(DF$AI - (DF$AI - tmpDFbt[, 2]) / sqrt(length(vectReps)) * CC,
  #                           function(lb){ max(0, lb) })
  # DF$BT_CIright_CC <- sapply(DF$AI + (tmpDFbt[, 3] - DF$AI) / sqrt(length(vectReps)) * CC,
  #                            function(ub){ min(1, ub) })

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


PerformBinTestAIAnalysisForConditionNPoint <- function(inDF, vectReps, pt = 0.5, binNObs=40, Q=0.95, EPS=1.3,
                                                       thr=NA, thrUP=NA, thrType="each", minDifference=NA){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for condition + point estimate to compare
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param pt A point to compare with
  #' @param binNObs Threshold on number of observations per bin
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param EPS An optional parameter to set a log window for coverage binning
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @param minDifference if specified, one additional column DAE is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to having non-overlapping CIs)
  #' @return A table of gene names, AIs + CIs, classification into genes demonstrating differential from point estimate AI and those that don't
  #' @examples
  #'

  fitDATA <- ComputeCorrConstantsForAllPairsReps(inDF, vectReps=vectReps, binNObs=binNObs,
                                                 EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType)

  vectRepsCombsCC <- sapply(fitDATA, function(fd){
    fd$fittedCC
  })

  RES <- PerformBinTestAIAnalysisForConditionNPoint_knownCC(inDF, vectReps=vectReps,
                                                            vectRepsCombsCC=vectRepsCombsCC,
                                                            pt=pt, binNObs=binNObs, Q=Q,
                                                            thr=thr, thrUP=thrUP, thrType=thrType,
                                                            minDifference=minDifference)

  return(list(CC = vectRepsCombsCC,
              FitDATA = fitDATA,
              Output = RES))
}


# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PERFORM DIFF CI(AI) ANALYSIS -- CONDITION & POINTS VECTOR
# ---------------------------------------------------------------------------------------

PerformBinTestAIAnalysisForConditionNPointVect_knownCC <- function(inDF, vectReps, vectRepsCombsCC, ptVect,
                                                                   binNObs=40, Q=0.95,
                                                                   thr=NA, thrUP=NA, thrType="each", minDifference=NA){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for condition + point estimate to compare
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param vectRepsCombsCC A vector of pairwise-computed correction constants for given replicates
  #' @param pt A point to compare with
  #' @param binNObs Threshold on number of observations per bin
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @param minDifference if specified, one additional column DAE is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to having non-overlapping CIs)
  #' @return A table of gene names, AIs + CIs, classification into genes demonstrating differential from point estimate AI and those that don't
  #' @examples
  #'

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

  # Bin test:
  tmpDFbt <- t(sapply(1:nrow(DF), function(i){
    if(is.na(matCOV[i]) | is.na(sumCOV[i]) | is.na(ptVect[i]) | sumCOV[i]==0) { return(c(NA,NA,NA)) }
    #print(paste(matCOV[i], sumCOV[i], ptVect[i]))
    BT <- binom.test(matCOV[i], sumCOV[i], alternative="two.sided", conf.level = Qbf, p=ptVect[i])
    c(BT$p.value, BT$conf.int[1], BT$conf.int[2])
  }))

  # DF$BT_pval = tmpDFbt[, 1]
  # DF$BT_CIleft = DF$AI - (DF$AI - tmpDFbt[, 2]) / sqrt(length(vectReps))
  # DF$BT_CIright = DF$AI + (tmpDFbt[, 3] - DF$AI) / sqrt(length(vectReps))
  #
  # DF$BT_CIleft_CC <- sapply(DF$AI - (DF$AI - tmpDFbt[, 2]) / sqrt(length(vectReps)) * CC,
  #                           function(lb){ max(0, lb) })
  # DF$BT_CIright_CC <- sapply(DF$AI + (tmpDFbt[, 3] - DF$AI) / sqrt(length(vectReps)) * CC,
  #                            function(ub){ min(1, ub) })

  DF$BT_pval = tmpDFbt[, 1]
  DF$BT_CIleft = tmpDFbt[, 2]
  DF$BT_CIright = tmpDFbt[, 3]

  DF$BT_CIleft_CC <- sapply(DF$AI - (DF$AI - tmpDFbt[, 2]) * CC,
                            function(lb){ max(0, lb) })
  DF$BT_CIright_CC <- sapply(DF$AI + (tmpDFbt[, 3] - DF$AI) * CC,
                             function(ub){ min(1, ub) })

  #DF$BT <- (DF$BT_pval < (1-Q)/nrow(na.omit(DF)))
  DF$BT <- sapply(1:nrow(DF), function(i){!(DF$BT_CIleft[i] <= ptVect[i] & DF$BT_CIright[i] >= ptVect[i])})
  # Find intersecting intervals > call them FALSE (non-rejected H_0)
  DF$BT_CC <- sapply(1:nrow(DF), function(i){!(DF$BT_CIleft_CC[i] <= ptVect[i] & DF$BT_CIright_CC[i] >= ptVect[i])})

  if (!is.na(minDifference))
  {
    DF$BT_CC_thrDiff <- (DF$BT_CC & (abs(DF$AI - ptVect[i]) >= minDifference))
  }

  return(cbind(DF, ptVect))
}


PerformBinTestAIAnalysisForConditionNPointVect <- function(inDF, vectReps, ptVect, binNObs=40, Q=0.95, EPS=1.3,
                                                           thr=NA, thrUP=NA, thrType="each", minDifference=NA){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for condition + point estimate to compare
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param pt A point to compare with
  #' @param binNObs Threshold on number of observations per bin
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param EPS An optional parameter to set a log window for coverage binning
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @param minDifference if specified, one additional column DAE is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to having non-overlapping CIs)
  #' @return A table of gene names, AIs + CIs, classification into genes demonstrating differential from point estimate AI and those that don't
  #' @examples
  #'

  fitDATA <- ComputeCorrConstantsForAllPairsReps(inDF, vectReps=vectReps, binNObs=binNObs,
                                                 EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType)

  vectRepsCombsCC <- sapply(fitDATA, function(fd){
    fd$fittedCC
  })

  RES <- PerformBinTestAIAnalysisForConditionNPointVect_knownCC(inDF, vectReps=vectReps,
                                                                vectRepsCombsCC=vectRepsCombsCC,
                                                                ptVect=ptVect, binNObs=binNObs, Q=Q,
                                                                thr=thr, thrUP=thrUP, thrType=thrType,
                                                                minDifference=minDifference)

  return(list(CC = vectRepsCombsCC,
              FitDATA = fitDATA,
              Output = RES))
}


# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PERFORM DIFF CI(AI) ANALYSIS -- 2 CONDITIONS
# ---------------------------------------------------------------------------------------

ComputeAICIs <- function(inDF, vectReps, vectRepsCombsCC,
                         Q=0.95, BF=T,
                         thr=NA, thrUP=NA, thrType="each"){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for condition + point estimate to compare
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param vectRepsCombsCC A vector of pairwise-computed correction constants for given replicates
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param BT Bonferroni correction, set False if Q is alredy corrected
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return A table of gene names, AIs + CIs, classification into genes demonstrating differential from point estimate AI and those that don't
  #' @examples
  #'

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
    BT <- binom.test(matCOV[i], sumCOV[i], alternative="two.sided", conf.level = Qbf)
    c(BT$conf.int[1], BT$conf.int[2])
  }))

  # DF$BT_CIleft = DF$AI - (DF$AI - tmpDFbt[, 1]) / sqrt(length(vectReps))
  # DF$BT_CIright = DF$AI + (tmpDFbt[, 2] - DF$AI) / sqrt(length(vectReps))
  #
  # DF$BT_CIleft_CC <- sapply(DF$AI - (DF$AI - tmpDFbt[, 1]) / sqrt(length(vectReps)) * CC,
  #                           function(lb){ max(0, lb) })
  # DF$BT_CIright_CC <- sapply(DF$AI + (tmpDFbt[, 2] - DF$AI) / sqrt(length(vectReps)) * CC,
  #                            function(ub){ min(1, ub) })

  DF$BT_CIleft = tmpDFbt[, 1]
  DF$BT_CIright = tmpDFbt[, 2]

  DF$BT_CIleft_CC <- sapply(DF$AI - (DF$AI - tmpDFbt[, 1]) * CC,
                            function(lb){ max(0, lb) })
  DF$BT_CIright_CC <- sapply(DF$AI + (tmpDFbt[, 2] - DF$AI) * CC,
                             function(ub){ min(1, ub) })
  return(DF)
}


PerformBinTestAIAnalysisForTwoConditions_knownCC <- function(inDF, vect1CondReps, vect2CondReps,
                                                             vect1CondRepsCombsCC, vect2CondRepsCombsCC,
                                                             Q=0.95,
                                                             thr=NA, thrUP=NA, thrType="each", minDifference=NA){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for condition + point estimate to compare
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vect1CondReps A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
  #' @param vect2CondReps A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
  #' @param vect1CondRepsCombsCC A vector of pairwise-computed correction constants for first condition's tech reps
  #' @param vect2CondRepsCombsCC A vector of pairwise-computed correction constants for second condition's tech reps
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @param minDifference if specified, one additional column DAE is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to having non-overlapping CIs)
  #' @return A table of gene names, AIs + CIs, classification into genes demonstrating differential from point estimate AI and those that don't
  #' @examples
  #'

  vectReps <- list(vect1CondReps, vect2CondReps)
  vectRepsCombsCC <- list(vect1CondRepsCombsCC,
                          vect2CondRepsCombsCC)

  DF_CI_divided <- lapply(1:2, function(i){
    ComputeAICIs(inDF, vectReps=vectReps[[i]], vectRepsCombsCC=vectRepsCombsCC[[i]],
                 Q=Q, BF=T, thr=thr, thrUP=thrUP, thrType=thrType)
  })

  DF <- merge(DF_CI_divided[[1]], DF_CI_divided[[2]], by = "ID")
  names(DF) <- c("ID",
                 paste0(names(DF_CI_divided[[1]])[-1], '_1'),
                 paste0(names(DF_CI_divided[[2]])[-1], '_2'))

  #################
  # p_bar = (DF$matCOV_1 + DF$matCOV_2)/(DF$sumCOV_1 + DF$sumCOV_2)
  # DF$Z <- (DF$AI_1 - DF$AI_2) / sqrt(p_bar * (1 - p_bar) * (1/DF$sumCOV_1 + 1/DF$sumCOV_2))


  # Find intersecting intervals > call them FALSE (non-rejected H_0)
  Qbft <- 1 - (1-Q)/nrow(na.omit(DF[, c("sumCOV_1", "matCOV_1", "sumCOV_2", "matCOV_2")]))
  DF_CI_divided_BFor2 <- lapply(1:2, function(i){
    ComputeAICIs(inDF, vectReps=vectReps[[i]], vectRepsCombsCC=vectRepsCombsCC[[i]],
                 Q=Qbft, BF=F,
                 thr=thr, thrUP=thrUP, thrType=thrType)
  })
  DF_BFor2 <- merge(DF_CI_divided_BFor2[[1]][, c("ID", "BT_CIleft", "BT_CIright")], DF_CI_divided_BFor2[[2]][, c("ID", "BT_CIleft", "BT_CIright")], by = "ID")
  names(DF_BFor2)[-1] <- c("BT_dCIleft_1", "BT_dCIright_1", "BT_dCIleft_2", "BT_dCIright_2")
  DF_BFor2_CC <- merge(DF_CI_divided_BFor2[[1]][, c("ID", "BT_CIleft_CC", "BT_CIright_CC")], DF_CI_divided_BFor2[[2]][, c("ID", "BT_CIleft_CC", "BT_CIright_CC")], by = "ID")
  names(DF_BFor2_CC)[-1] <- c("BT_dCIleft_CC_1", "BT_dCIright_CC_1", "BT_dCIleft_CC_2", "BT_dCIright_CC_2")

  #print(DF)

  DF <- merge(DF, DF_BFor2, by = "ID")
  DF$BT <- (DF$BT_dCIright_2 < DF$BT_dCIleft_1 | DF$BT_dCIright_1 < DF$BT_dCIleft_2)

  DF <- merge(DF, DF_BFor2_CC, by = "ID")
  DF$BT_CC <- (DF$BT_dCIright_CC_2 < DF$BT_dCIleft_CC_1 | DF$BT_dCIright_CC_1 < DF$BT_dCIleft_CC_2)

  #print(DF)

  if (!is.na(minDifference))
  {
    DF$BT_CC_thrDiff <- (DF$BT_CC & (abs(DF$AI_1 - DF$AI_2) >= minDifference))
  }

  return(DF)
}


PerformBinTestAIAnalysisForTwoConditions <- function(inDF, vect1CondReps, vect2CondReps, binNObs=40, Q=0.95, EPS=1.3,
                                                     thr=NA, thrUP=NA, thrType="each", minDifference=NA){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for condition + point estimate to compare
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vect1CondReps A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
  #' @param vect2CondReps A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
  #' @param binNObs Threshold on number of observations per bin
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param EPS An optional parameter to set a log window for coverage binning
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @param minDifference if specified, one additional column DAE is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to having non-overlapping CIs)
  #' @return A table of gene names, AIs + CIs, classification into genes demonstrating differential from point estimate AI and those that don't
  #' @examples
  #'

  fitDATA1Cond <- ComputeCorrConstantsForAllPairsReps(inDF, vectReps=vect1CondReps, binNObs=binNObs,
                                                      EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType)
  fitDATA2Cond <- ComputeCorrConstantsForAllPairsReps(inDF, vectReps=vect2CondReps, binNObs=binNObs,
                                                      EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType)

  vect1CondRepsCombsCC <- sapply(fitDATA1Cond, function(fd){
    fd$fittedCC
  })
  vect2CondRepsCombsCC <- sapply(fitDATA2Cond, function(fd){
    fd$fittedCC
  })
  print(paste(vect1CondRepsCombsCC, vect2CondRepsCombsCC))

  RES <- PerformBinTestAIAnalysisForTwoConditions_knownCC(inDF, vect1CondReps=vect1CondReps, vect2CondReps=vect2CondReps,
                                                          vect1CondRepsCombsCC=vect1CondRepsCombsCC,
                                                          vect2CondRepsCombsCC=vect2CondRepsCombsCC,
                                                          Q=Q,
                                                          thr=thr, thrUP=thrUP, thrType=thrType,
                                                          minDifference=minDifference)

  return(list(CC = list(vect1CondRepsCombsCC, vect2CondRepsCombsCC),
              FitDATA = list(fitDATA1Cond, fitDATA2Cond),
              Output = RES))
}

nameColumns <- function(exp_n, rep_n)  {
  #' Helper function to quickly rename columns in geneCountTab dataframe
  #'
  #' @param exp_n Experiment number
  #' @param rep_n Number of replicates for the experiment
  #' @return Vector with names
  #' @examples colnames(geneCountTab)[2:13] <- c(nameColumns(1,2), nameColumns(2,2), nameColumns(3,2))
  #'

  paste0("exp", rep(exp_n, 2*rep_n), "_rep", rep(1:rep_n, each = 2), "_", rep(c("ref", "alt"), rep_n))
}
