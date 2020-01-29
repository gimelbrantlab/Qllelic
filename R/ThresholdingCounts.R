#' Takes table with gene names and counts and returns table, where all genes that under threshold have NA coverage. Can be restricted to particular reps.
#'
#' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns
#' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df)
#' @param thr An optional parameter for a threshold on gene coverage (default = NA)
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#'
#' @return Table with masked with NA undercovered genes
#' @export
#'
ThresholdingCounts <- function(df, reps=NA, thr=NA, thrUP=NA, thrType="each"){

  options(stringsAsFactors = FALSE)
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
