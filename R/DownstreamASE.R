# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PROCESS TABLES WITH COUNTs/AIs, FIND MAE GENES
# ---------------------------------------------------------------------------------------

DiffAIplusFig <- function(inDF, CondReps_expA, CondReps_expB,
                          CondReps_expA_name, CondReps_expB_name,
                          expA_CC=NA, expB_CC=NA,
                          thr_coverage=NA, minDifference=0.1, lm=T, save=F) {
  #' A function that gets a table with ref & alt counts per gene/SNP for each replicate
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param CondReps_expA A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
  #' @param CondReps_expB A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
  #' @param CondReps_expA_name Name of the first condition
  #' @param CondReps_expB_name Name of the second condition
  #' @param thr_coverage optional parameter, specifies coverage threshold
  #' @param minDifference minimal difference between AIs in 2 conditions to call them differentially allelically expressed
  #' @param save If true, save the results to table
  #' @return inDF table with 2 additional columns classifying genes into differentially allelically expressed (DAE) and not + figure
  #' @examples
  #'
  if ((!is.na(expA_CC))&(!is.na(expB_CC))) {
    df_tmp <- PerformBinTestAIAnalysisForTwoConditions_knownCC(inDF,
                                                               vect1CondReps = CondReps_expA,
                                                               vect2CondReps = CondReps_expB,
                                                               vect1CondRepsCombsCC = expA_CC,
                                                               vect2CondRepsCombsCC = expB_CC,
                                                               Q=0.95,
                                                               thr=thr_coverage,
                                                               minDifference=minDifference)
  }
  else {
    df_tmp <- PerformBinTestAIAnalysisForTwoConditions(inDF,
                                                       vect1CondReps = CondReps_expA,
                                                       vect2CondReps = CondReps_expB,
                                                       Q=0.95,
                                                       thr=thr_coverage,
                                                       minDifference=minDifference)
  }
  fig_compare_tmp <- ggplot(df_tmp, aes(x = AI_1, y = AI_2, col = BT_CC_thrDiff)) +
    geom_point(size=0.5) +
    theme_bw() +
    xlab(paste0("AI, ",CondReps_expA_name)) +
    ylab(paste0("AI, ",CondReps_expB_name)) +
    scale_color_manual(name="Differential AI", labels=c("FALSE", "TRUE"), values=c("gray", "red")) +
    coord_fixed() +
    theme(legend.position = "None")
  if (lm)
    fig_compare_tmp <- fig_compare_tmp + geom_smooth(method="lm")
  if (save) {
    AI_table_with_CIs <- df_tmp[,c(1,4,21,22,11,23,24,25,26)]
    colnames(AI_table_with_CIs) <- c("ensembl_gene_id",
                                     paste0("AI_",CondReps_expA_name),
                                     paste0("AI_CI_left_",CondReps_expA_name),
                                     paste0("AI_CI_right_",CondReps_expA_name),
                                     paste0("AI_",CondReps_expB_name),
                                     paste0("AI_CI_left_",CondReps_expB_name),
                                     paste0("AI_CI_right_",CondReps_expB_name),
                                     "CI_diff",
                                     "CI_plus_minDiff")
    write.table(AI_table_with_CIs,
                file=paste0("~/Dropbox (Partners HealthCare)/MAE screen/Drugs/Manuscript/Tables_SV/add_clones_5aza_treated/AI_", CondReps_expB_name,"_vs_",CondReps_expA_name,"_min_diff_",minDifference,"_cov_threshold_",thr_coverage,"_upd_May.txt"), quote = F, row.names = F)
  }
  return(list(df_tmp, fig_compare_tmp))
}



findMAE <- function(x) {
  #' A function that returns classification based on AI
  #'
  #' @param x Vector with ASE classification per clone
  #' @return overall MAE status among clones (by Gendrel classification)
  #' @examples
  #'
  mon <- 0
  if (length(x)==1) {
    mon=NA
  }
  else {
    not_nm_count <- length(x) - sum(x=="nd")
    if (not_nm_count==0) mon="nd"
    else if ((sum(x=="CAST_monoallelic")+sum(x=="CAST_biased")==not_nm_count)|(sum(x=="129_monoallelic")+sum(x=="129_biased")==not_nm_count)) mon="gen_sk"
    else if ((sum(x=="CAST_monoallelic")>0)|(sum(x=="129_monoallelic")>0)) mon="monoallelic"
    else if ((sum(x=="CAST_biased")>0)|(sum(x=="129_biased")>0)) mon="biased"
    else if (sum(x=="biallelic")==not_nm_count) mon="biallelic"
    else mon="other"
  }
  return(mon)
}

isMAE_test_CI <- function(x) {
  #' A function that tests for allelic imbalance
  #'
  #' @param x Vector two entries: first is the output of CIs test (T if genes are differentially monoallelically expressed), second is absolute AI value
  #' @return overall MAE status among clones (by Gendrel classification)
  #' @examples
  #'
  thr <- x[3]
  if ((is.na(x[1]))|(is.na(x[2]))) {
    return("nd")
  }
  else {
    if ((x[1])&((x[2]>=thr))) {
      return("129_monoallelic")
    }
    else if ((x[1])&((x[2]<=(1-thr)))) {
      return("CAST_monoallelic")
    }
    else if ((x[1])&((x[2]>=0.5))) {
      return("129_biased")
    }
    else if ((x[1])&((x[2]<0.5))) {
      return("CAST_biased")
    }
    else {return("biallelic")}
  }
}
