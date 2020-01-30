#' MeanCoverage
#'
#' Calculates mean coverage mat+pat among given replicates.
#'
#' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
#' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df).
#' @param thr An optional parameter for a threshold on gene coverage (default = NA)
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#'
#' @return Df with names and Mean coverage mat+pat among given replicates.
#' @export
#'
#' @examples MeanCoverage(allelicCountsTable, reps=c(3,4), thr=8)
#'
MeanCoverage <- function(df, reps=NA, thr=NA, thrUP=NA, thrType="each"){
  ddf <- ThresholdingCounts(df, reps, thr, thrUP, thrType)

  res_df <- data.frame(df[, 1], meanCOV = rowMeans(ddf[, -1])*2, stringsAsFactors = FALSE)
  names(res_df)[1] <- names(df)[1]

  return(res_df)
}
