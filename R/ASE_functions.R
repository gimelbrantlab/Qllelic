# _______________________________________________________________________________________

options(stringsAsFactors = FALSE)
# _______________________________________________________________________________________

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: READ FILES AND CREATE UNIFORM TABLES
# ---------------------------------------------------------------------------------------
#                 TODO: 1. Kallisto functions
# ---------------------------------------------------------------------------------------

GetGatkPipelineTabs <- function(inFiles, nReps, multiple = TRUE, chrom = F){
  #' (GATK pipeline) Concatenate vertically (uniting merge) tables from inFiles.
  #'
  #' @param inFiles A vector of full pathes to files with gene names and allelic counts
  #' @param nReps A vector of numbers, each entry is a number of replicates in the corresponding file
  #' @param multiple Parameter defining if multiple input files are used, default set to TRUE
  #' @param chrom Parameter defining if the resulting table includes chromosome column, default set to FALSE
  #' @return A concatenated table with means/counts, each row corresponds to a gene
  #' @examples
  #'
  # TODO : change naming here
  if (multiple) {
    df <- Reduce(function(x,y){merge(x,y,by="ensembl_gene_id")},
                 lapply(1:length(inFiles), function(x){
                   df0 <- read_delim(inFiles[x], delim="\t", escape_double = FALSE, trim_ws = TRUE)
                   return(df0[,c(1:(2*nReps[x]+1))])
                 }
                 )
    )
    # TODO add chrom option for multiple files
  }
  else {
    df <- read_delim(inFiles, delim="\t", escape_double = FALSE, trim_ws = TRUE)
    df_chrom <- df[,c("ensembl_gene_id","chr")]
    df <- df[,c(1:(2*sum(nReps)+1))]
    df <- as.data.frame(df)
    if (chrom) {
      df <- merge(df, df_chrom, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
    }
  }
  nameColumns <- function(rep_n)  {
    paste0("rep", rep(1:rep_n, each = 2), rep(c("_ref", "_alt"), rep_n))
  }
  if (multiple) {
    names(df) <- c("ensembl_gene_id",
                   paste0("rep", rep(1:sum(nReps), each=2), c("_ref", "_alt")))
    if (chrom) {
      names(df)[ncol(df)] <- "chr"
    }
  }
  else {
    names(df) <- c("ensembl_gene_id", unlist(sapply(nReps, nameColumns)))
    if (chrom) {
      names(df)[ncol(df)] <- "chr"
    }
  }
  return(df)
}

GetGatkPipelineSNPTabs <- function(inFiles, nReps){
  #' (GATK pipeline) Concatenate vertically (uniting merge) tables from inFiles
  #'
  #' @param inFiles A vector of pathes to files
  #' @param nReps A vector of numbers, each -- number of replicates in corresponding file
  #' @return A technical replicates-concatenated table with SNPs, rows corresponds to SNPs
  #' @examples
  #'
  df <- Reduce(function(x,y){merge(x,y,by="ensembl_gene_id")},
               lapply(1:length(inFiles), function(x){
                 df <- read_delim(inFiles[x], delim="\t", escape_double = FALSE, trim_ws = TRUE)
                 cs <- which(sapply(names(df), function(x){
                   (grepl("rep", x, fixed=TRUE) & grepl("aggr", x, fixed=TRUE))
                 }
                 )
                 )
                 snpdf <- apply(df[, cs], 2, function(c){as.numeric(unlist(strsplit(c, ', ')))})
                 snpdf <- data.frame(ensembl_gene_id = paste0("SNP", 1:nrow(snptab_list[[1]])),
                                     snpdf)
                 return(snpdf)
               }
               )
  )
  names(df) <- c("ensembl_gene_id",
                 paste0("rep", rep(1:sum(nReps), each=2), c("_ref", "_alt")))
  return(df)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: ALLELIC IMBALANSE AND MEAN COVERAGE
# ---------------------------------------------------------------------------------------

ThresholdingCounts <- function(df, reps=NA, thr=NA, thrUP=NA, thrType="each"){
  #' Takes table with gene names and counts and returns table, where all genes that under threshold have NA coverage. Can be restricted to particular reps.
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df)
  #' @param thr An optional parameter for a threshold on gene coverage (default = NA)
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return Table with masked with NA undercovered genes
  #' @examples
  #'

  # Taking replicates:
  if(anyNA(reps)){
    reps <- 1:(ncol(df)%/%2)
    ddf  <- df
  } else {
    cs  <- sort(c(sapply(reps, function(x){c(x*2, x*2+1)}))) # columns numbers
    ddf <- df[, c(1, cs)]
  }

  # Thresholding:
  if(!anyNA(thr)){
    if(thrType == "each"){
      # Masking with NA gene coverage info replicate-specific if it is < thr
      greaterThanThr <- (sapply(1:length(reps), function(x){
        rep_thr_info <- rowSums(ddf[, c(2*x, 2*x+1)])
        cbind(rep_thr_info, rep_thr_info)
      }) >= thr)
      greaterThanThr[greaterThanThr==FALSE] <- NA
      ddf[, -1] <- ddf[, -1] * greaterThanThr
    } else if(thrType == "average"){
      # Masking with NA entire gene coverage info for all replicates if average coverage < thr:
      greaterThanThr <- (rowMeans(ddf[, -1])*2 >= thr)
      greaterThanThr[greaterThanThr==FALSE] <- NA
      ddf[, -1] <- ddf[, -1] * greaterThanThr
    }
  }

  #Same for UPPER Threshold:
  if(!anyNA(thrUP)){
    if(thrType == "each"){
      # Masking with NA gene coverage info replicate-specific if it is > thrUP
      lesserThanThrUP <- (sapply(1:length(reps), function(x){
        rep_thr_info <- rowSums(ddf[, c(2*x, 2*x+1)])
        cbind(rep_thr_info, rep_thr_info)
      }) <= thrUP)
      lesserThanThrUP[lesserThanThrUP==FALSE] <- NA
      ddf[, -1] <- ddf[, -1] * lesserThanThrUP
    } else if(thrType == "average"){
      # Masking with NA entire gene coverage info for all replicates if average coverage > thr:
      lesserThanThrUP <- (rowMeans(ddf[, -1])*2 <= thrUP)
      lesserThanThrUP[lesserThanThrUP==FALSE] <- NA
      ddf[, -1] <- ddf[, -1] * lesserThanThrUP
    }
  }

  return(ddf)
}

CountsToAI <- function(df, reps=NA, meth="meanOfProportions", thr=NA, thrUP=NA, thrType="each"){
  #' Calculates allelic imbalances from merged counts over given replicates (ai(sum_reps(gene)))
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df)
  #' @param meth An optional parameter for method to use, either sum(m)/sum(p), or sum(m/p) (default = sum(m)/sum(p))
  #' @param thr An optional parameter for a threshold on gene coverage (default = NA)
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return Df with names and mean(mean(m_1,...,m_6))_SNP / mean(mean(m_1+p_1,...,m_6+p6))_SNP
  #' @examples
  #'

  ddf <- ThresholdingCounts(df, reps, thr, thrUP, thrType)

  if(ncol(ddf) == 3){ # if 1 replicate
    p <- (ddf[, 2]/rowSums(ddf[, -1]))
  } else {             # if more than 1 replicates
    if (meth == "mergedToProportion") {
      ref <- rowSums(ddf[, seq(2, ncol(ddf), 2)])
      all <- rowSums(ddf[, 2:ncol(ddf)])
      p   <- ref/all
    } else if(meth == "meanOfProportions"){
      aitab <- sapply(1:(ncol(ddf)%/%2), function(i){
        ddf[, i*2]/(ddf[, i*2]+ddf[, i*2+1])
      })
      p <- rowMeans(aitab)
    }
  }
  p[is.nan(p)] <- NA

  res_df <- data.frame(df[, 1], AI = p)
  names(res_df)[1] <- names(df)[1]

  return(res_df)
}

MeanCoverage <- function(df, reps=NA, thr=NA, thrUP=NA, thrType="each"){
  #' Calculates mean coverage mat+pat among given replicates.
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @param thr An optional parameter for a threshold on gene coverage (default = NA)
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return Df with names and Mean coverage mat+pat among given replicates.
  #' @examples
  #'

  ddf <- ThresholdingCounts(df, reps, thr, thrUP, thrType)

  res_df <- data.frame(df[, 1], meanCOV = rowMeans(ddf[, -1])*2)
  names(res_df)[1] <- names(df)[1]

  return(res_df)
}

MergeSumCounts <- function(df, reps=NA, thr=NA, thrUP=NA, thrType="each"){
  #' Creates a table of sums of maternal and paternal count for given replicates.
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @param thr An optional parameter for a threshold on gene coverage (default = NA)
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return Table with names and sums of mat and pat coverages among given replicates.
  #' @examples
  #'

  ddf <- ThresholdingCounts(df, reps, thr, thrUP, thrType)

  if(ncol(ddf) == 3){
    res_df <- data.frame(ddf[, 1], ref_reps = ddf[, 2], alt_reps = ddf[, 3])
  } else {
    res_df <- data.frame(ddf[, 1],
                         ref_reps = rowSums(ddf[, seq(2, ncol(ddf), 2)]),
                         alt_reps = rowSums(ddf[, seq(3, ncol(ddf), 2)]))
  }
  names(res_df)[1] <- names(ddf)[1]

  return(res_df)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: NAMES REWRITING
# ---------------------------------------------------------------------------------------

NumToDoulbledigitChar <- function(x){
  #' Creates double-digit character from a number or a vector of numbers (up t0 99)
  #'
  #' @param x A number from 0 to 99
  #' @return Double-digit character (or vector) from a number ('x' or '0x')
  #' @examples
  #'
  if(length(x)==1){ x = c(x) }
  doubledigitsChars <- sapply(1:length(x), function(i){
    if(x[i]<=9) { return(paste0('0',x[i])) }
    else { return(as.character(x[i])) }
  })
  return(doubledigitsChars)
}

BuildDesign <- function(experimentNames, techReps, corrConst=NA){
  #' Creates a design matrix for the experiment
  #'
  #' @param experimentNames Vector with names of the experiments
  #' @param techReps Vector with number of technical replicates in each experiment
  #' @param corrConst Optional, Vector with correction constants for each experiment
  #' @return Dataframe with experiments numbered and numbers of columns
  #' @examples
  #'
  rowsSp <- data.frame(matrix(lapply(1:length(techReps), function(x){(2*sum(techReps[1:x-1])+1):(2*sum(techReps[1:x]))}), nrow=length(techReps), byrow=T),stringsAsFactors=FALSE)
  colnames(rowsSp) <- "replicateCols"
  colExp <- data.frame(matrix(lapply(1:length(techReps), function(x){(sum(techReps[1:x])-techReps[x]+1):(sum(techReps[1:x]))}), nrow=length(techReps), byrow=T),stringsAsFactors=FALSE)
  colnames(colExp) <- "replicateNums"
  designMatrix <- cbind(experimentNames, techReps, rowsSp, colExp)
  colnames(designMatrix) <- c("experimentNames", "techReps", "replicateCols", "replicateNums")
  if (!sum(is.na(corrConst))) {
    designMatrix <- cbind(designMatrix, corrConst)
  }
  return(designMatrix)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PAIRED COMPARISONS
# ---------------------------------------------------------------------------------------

CreateDeltaAIForAPairRepsDF <- function(df, reps, thr=NA, thrUP=NA, thrType="each"){
  #' Creates a tab with gene mean coverage and AI deltas for a pair of technical replicates
  #'
  #' @param tab A dataframe of genes/transcripts and parental counts for technical replicates in columns
  #' @param reps A parameter for setting 2 replicates for consideration
  #' @param thr An optional parameter for a threshold on gene coverage (default = NA)
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return A tab with names, gene mean coverage and AI deltas for a pair of technical replicates
  #' @examples
  #'
  ai1 <- CountsToAI(df, reps[1], thr=thr, thrUP=thrUP, thrType=thrType)$AI
  ai2 <- CountsToAI(df, reps[2], thr=thr, thrUP=thrUP, thrType=thrType)$AI
  ddf <- data.frame(df[, 1],
                    deltaAI = ai1 - ai2,
                    AI1 = ai1, AI2 = ai2,
                    MeanCov = MeanCoverage(df, reps, thr=thr, thrUP=thrUP, thrType=thrType)$meanCOV)
  names(ddf)[1] <- names(df)[1]

  return(na.omit(ddf))
}

CreateMergedDeltaAIPairwiseDF <- function(df, reps=NA, what="noname",
                                          thr=NA, thrUP=NA, thrType="each"){
  #' Creates a table of parwise comparisons for technical replicates in given table
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df)
  #' @param what A name (default = "noname")
  #' @param thr An optional parameter for a threshold on gene coverage (default = NA)
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return A table of parwise comparisons for technical replicates in given table
  #' @examples
  #'

  # Take all combinations from given replicates:
  if(anyNA(reps)) {
    reps <- 1:(ncol(df)%/%2)
  }
  combs <- combn(reps, 2, function(x){x})

  # Perform paiwise ai comparison for all replicate combinations:
  ddfPairs <- do.call(rbind, lapply(1:ncol(combs), function(y){
    ddf          <- CreateDeltaAIForAPairRepsDF(df, combs[,y], thr, thrUP, thrType)
    ddf$deltaAI  <- abs(ddf$deltaAI)
    ddf$what     <- what
    ddf$i   <- combs[1,y]
    ddf$j   <- combs[2,y]
    ddf$ij  <- paste(NumToDoulbledigitChar(combs[1,y]),
                     'vs',
                     NumToDoulbledigitChar(combs[2,y]))
    return(ddf)
  }))

  return(ddfPairs)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: BINNING AND QUARTILLING
# ---------------------------------------------------------------------------------------

CreateObservedQuantilesDF <- function(df, P, ep, logbase=T, coverageLimit, group=''){
  #' Creates a table with quantiles and numbers of bins for technical replicates for a given table, binned into log intervals
  #'
  #' @param df A dataframe - output of CreateMergedDeltaAIPairwiseDF()
  #' @param P A vector of %-quartiles
  #' @param logbase The binary parameter if we deal with log base (default = T)
  #' @param ep The log base for binning (0^b, 1^b, ...)
  #' @param coverageLimit Gene coverage limit for consideration
  #' @param group An optional name (default = '')
  #' @return A table with quantiles and numbers of bins for technical replicates for a given table
  #' @examples
  #'
  if(logbase){
    covIntervalsStarts <- unique(floor(ep**(0:log(coverageLimit, base=ep))))
  } else {
    covIntervalsStarts <- seq(0, coverageLimit-ep, ep)
  }
  lenCIS       <- length(covIntervalsStarts)
  covIntervals <- c(covIntervalsStarts, coverageLimit)


  ddf <- do.call(rbind, lapply(P, function(p){ # [for all quartile (%)]:
    # [for all coverage bins]:
    ddfP <- data.frame(coverageBin = covIntervalsStarts,
                       deltaAI = sapply(1:lenCIS, function(i){
                         dai <- df[df$MeanCov >= covIntervals[i] &
                                     df$MeanCov <  covIntervals[i+1], ]$deltaAI
                         quantile(dai, p, na.rm = T)
                       }),
                       binNObservations = sapply(1:lenCIS, function(i){
                         nrow(df[df$MeanCov >= covIntervals[i] &
                                   df$MeanCov <  covIntervals[i+1], ])
                       }),
                       Q = p)
    return(ddfP)
  })
  )
  ddf$group <- as.factor(paste(group, ddf$Q))
  ddf$Q     <- as.factor(ddf$Q)
  return(ddf)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: FITTING LM
# ---------------------------------------------------------------------------------------

# Was: fit_lm_intercept_how_we_want_morethan(inDF, N_obs_bin, morethan=10)
#
# TODO: ADD LESSTEN?
#
FitLmIntercept <- function(inDF, binNObs, morethan = 10, logoutput = TRUE){
  #' Fits linear model to logarithmic data and outputs intercept for the model with slope=1/2 restriction
  #
  #' @param inDF A dataframe - output of CreateObservedQuantilesDF()
  #' @param N_obs_bin Threshold on number of observations per bin
  #' @param morethan Theshold on gene coverage for lm (default = 10)
  #' @param logoutput Return log intercept? (default = true)
  #' @return lm intercept or log2(lm intercept)
  #' @examples
  #'
  df <- inDF[, c('coverageBin','deltaAI','binNObservations')]
  df <- df[df$coverageBin > morethan &
             df$binNObservations > binNObs, ]
  df[, c('coverageBin','deltaAI')] <- log2(df[, c('coverageBin','deltaAI')])
  loglm <- lm(data = df, deltaAI ~ offset(-0.5*coverageBin), na.action=na.exclude)$coefficients
  if(logoutput){
    return(loglm[1])
  } else {
    return(2**loglm[1])
  }
}


# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: AI CI ESTIMATES
# ---------------------------------------------------------------------------------------
CreatePMforAI <- function(dfInt, dfAI, dfCov){
  #' Input: two data frames, one including replicate-gene coverages, one with constants for each technical replicates pair
  #'
  #' @param dfInt A table with column "linInt" of correction constants for each replicates combination (rows), the order of rows should be consistent with columns in dfCov, s.t. pairs are alphabetically ordered
  #' @param dfAI A table with columns for each technical replicate, the rows correspond to genes, the values are AI
  #' @param dfCov A table with columns for each technical replicate, the rows correspond to genes, the values are coverage
  #' @return Plus-minus intervals to determine AI Confidence Intervals for each gene
  #' @examples
  #'
  covSumsCombs <- combn(1:ncol(dfCov), 2, function(x){rowSums(dfCov[, x])})
  covSumsCombs[rowMeans(is.na(dfAI))>0, ] <- NA
  invertCovSumsCombs <- 1 / covSumsCombs
  if(ncol(dfCov) == 2){
    qres = apply(invertCovSumsCombs, 1, function(c){c * dfInt$linInt**2})
  } else if(ncol(dfCov) > 2){
    qres = rowSums(t(apply(invertCovSumsCombs, 1, function(c){c * dfInt$linInt**2})))
  }
  return(0.5/(ncol(dfCov)*(ncol(dfCov)-1))*sqrt(qres))
}

CreateCIforAI <- function(dfInt, dfCounts, condName="Condition", thr=NA, thrUP=NA, thrType="each"){
  #' Input: two data frames, one including replicate-gene mat|pat coverages, one with constants for each technical replicates pair
  #'
  #' @param dfInt A table with column "linInt" of correction constants for each replicates combination (rows), the order of rows should be consistent with columns in dfCounts, s.t. pairs are alphabetically ordered
  #' @param dfCounts A table with pairs of columns for each technical replicate, the rows correspond to genes, the values are mat|pat coverages
  #' @param condName An optional parameter; one-word name for condition
  #' @param thr An optional parameter for a threshold on gene coverage (default = NA)
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return 5-column df: ID, AI, Plus-minus intervals, AI Confidence Intervals left and right ends
  #' @examples
  #'
  dfAI  <- sapply(1:(ncol(dfCounts)%/%2), function(i){
    CountsToAI(dfCounts, reps=i, thr=thr, thrUP=thrUP, thrType=thrType)$AI
  })
  dfCov <- sapply(1:(ncol(dfCounts)%/%2), function(i){
    dfCounts[, (2*i)] + dfCounts[, (2*i+1)]
  })

  df <- data.frame(
    ID        = dfCounts[,1],
    condition = condName,
    meanCov   = rowMeans(dfCov),
    meanAI    = CountsToAI(dfCounts, meth="meanOfProportions", thr=thr, thrUP=thrUP, thrType=thrType)$AI,
    pm        = CreatePMforAI(dfInt, dfAI, dfCov)
  )
  df$meanAILow  <- sapply(df$meanAI-df$pm, function(x){max(0, x)})
  df$meanAIHigh <- sapply(df$meanAI+df$pm, function(x){min(1, x)})
  return(df)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PERFORM QUANTILE ANALYSIS
# ---------------------------------------------------------------------------------------
PerformDAIQuantilesAnalysis <- function(inDF, vectReps, condName="Condition",
                                        Q=0.95, EPS=1.3, thr=NA, thrUP=NA, thrType="each"){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param condName An optional parameter; one-word name for condition
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param EPS An optional parameter to set a log window for coverage binning
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @param fullOUT Set true if you want full output with all computationally-internal dfs
  #' @return A table of quantiles in coverage bins
  #' @examples
  #'

  # Take subtable:
  dfCondition <- inDF[, sort(c(1, vectReps*2, vectReps*2+1))]

  # Create pairvise AI differences for all techreps pairs:
  deltaAIPairwiseDF <- CreateMergedDeltaAIPairwiseDF(df=dfCondition, what=condName, thr=thr, thrUP=thrUP, thrType=thrType)
  deltaAIPairwiseDF$group <- paste(condName, deltaAIPairwiseDF$ij)

  # Count quantiles for Mean Coverage bins:
  observedQuantilesDF <- do.call(rbind,
                                 lapply(unique(deltaAIPairwiseDF$group),
                                        function(gr){
                                          df  <- deltaAIPairwiseDF[deltaAIPairwiseDF$group == gr, ]
                                          res <- CreateObservedQuantilesDF(df=df,
                                                                           P=Q, ep=EPS, logbase=T,
                                                                           coverageLimit=quantile(deltaAIPairwiseDF$MeanCov, 0.995),
                                                                           group=gr)
                                        }
                                 )
  )
  observedQuantilesDF$condition <- condName
  observedQuantilesDF$ij <- sapply(as.character(observedQuantilesDF$group),
                                   function(x){paste(unlist(strsplit(x, ' '))[2:4], collapse=' ')})

  return(observedQuantilesDF)
}

# .......................................................................................
# Aftercomments:
# _______________________________________________________________________________________


# THE END
