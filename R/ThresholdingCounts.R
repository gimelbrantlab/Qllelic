#' ThresholdingCounts
#'
#' Takes allelic counts table and returns table, where all genes that don't pass a given coverage threshold have NA coverage. Can be restricted to particular replicates.
#'
#' @param df Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)
#' @param reps Optional (default=NA, all replicates), a vector of replicate numbers for which the analysis should be applied
#' @param thr Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis
#' @param thrUP Optional (default=NA), threshold for max gene coverage (default = NA)
#' @param thrType Optional (default = "each", also can be "average" for average coverage on replicates), threshold type
#'
#' @return Allelic counts table with masked with NA undercovered genes, for selected replicates
#' @export
#'
#' @examples ThresholdingCounts(df = allelicCountsTable, reps = c(1,2), thr = 10)
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
