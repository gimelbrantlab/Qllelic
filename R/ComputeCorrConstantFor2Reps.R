#' ComputeCorrConstantFor2Reps
#'
#' Computes QCC for one pair of replicates.
#'
#' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
#' @param reps A vector of 2 of replicate numbers that should be considered
#' @param binNObs Threshold on number of observations per bin
#' @param fitCovThr Threshold on coverage for genes that will be included in Beta-Bin fitting
#' @param EPS An optional parameter to set a log window for coverage binning
#' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#'
#' @return List with (1) fitted QCC ($fittedCC) and (2) a table with proportions of observed to expected quantiles per coverage bin ($QObsExpPropsTable).
#' @export
#'
#' @importFrom stats "lm" "na.exclude" "na.omit" "quantile" "rbeta" "rbinom"
#' @importFrom utils "combn"
#'
ComputeCorrConstantFor2Reps <- function(inDF, reps, binNObs=40, fitCovThr=50,
                                        EPS=1.3, thr=NA, thrUP=NA, thrType="each"){

  options(stringsAsFactors = FALSE)
  ##-------------------------------------------------------------------------------------------------------------------------------------
  ## 1. Compute beta-binomial parameters for merged AI distribution per each coverage bin:
  ##-------------------------------------------------------------------------------------------------------------------------------------
  ##     1.1. Table with counts and ai:
  ##
  if(EPS <= 1){ EPS = 1.001 }
  meancov =  MeanCoverage(inDF, reps=reps, thr=thr, thrUP=thrUP, thrType=thrType)$meanCOV
  df_unit_info = data.frame(ID = inDF[, 1],
                            AI_mean = CountsToAI(inDF, reps=reps, meth="meanOfProportions", thr=thr, thrUP=thrUP, thrType=thrType)$AI,
                            AI_merged = CountsToAI(inDF, reps=reps, meth="mergedToProportion", thr=thr, thrUP=thrUP, thrType=thrType)$AI,
                            AI_1 = CountsToAI(inDF, reps=reps[1], thr=thr, thrUP=thrUP, thrType=thrType)$AI,
                            AI_2 = CountsToAI(inDF, reps=reps[2], thr=thr, thrUP=thrUP, thrType=thrType)$AI,
                            mCOV = meancov,
                            COV = meancov*2,
                            binCOV = ceiling(EPS**(ceiling(log(meancov, base=EPS)) - 1/2)))
  # ceiling(log(min(meancov[!is.na(meancov)]), base=EPS))
  binCOVs = sort(unique(
    ceiling(EPS**c(-Inf, 0:ceiling(log(max(meancov[!is.na(meancov)]), base=EPS)) - 1/2))
  ))
  df_covbinsnum = data.frame(binCOV = binCOVs,
                             binNUM = sapply(binCOVs, function(x){
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
                                       df_covbinsnum$binCOV >= max(fitCovThr, thr)]

  print(paste(length(covbinsGthr), "COVERAGE BINS"))

  df_betabin_res = lapply(1:length(covbinsGthr), function(i){
    coverage = covbinsGthr[i]
    df = df_unit_info[df_unit_info$binCOV == coverage, ]
    observations = round(df$AI_merged * coverage*2)

    fitdat = MixBetaBinomialFit(initials, coverage*2, observations)
    fitres = fitdat[[1]]
    dfres = data.frame(No = i,
                       coverage = coverage,
                       x_weight_predicted = fitres[1],
                       x_alpha_predicted = fitres[2],
                       y_alpha_predicted = fitres[3],
                       algm_steps = fitres[4])
    print(paste(i, "|", "meanCOV:", covbinsGthr[i], ",", "#STEPS:", dfres$algm_steps, sep="    "))
    return(list(dfres, fitdat[[2]]))
  })
  df_betabin_params = do.call(rbind, lapply(1:length(covbinsGthr), function(i){
    df_betabin_res[[i]][[1]]
  }))
  fittabinfo = lapply(1:length(covbinsGthr), function(i){
    df_betabin_res[[i]][[2]]
  })

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
    # infoTable = df_unit_info,
    # betaBinTable = df_betabin_params,
    # betaBinFitInfo = fittabinfo,
    # dAIBetaBinFit = lst_dai_betabinfit,
    # dAIObserved = lst_dai_observed,
    QObsExpPropsTable = df_observed_expected_quantile_proportions,
    fittedCC = fittedCC
  ))

}
