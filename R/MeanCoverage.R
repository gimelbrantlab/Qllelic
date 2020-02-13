#' MeanCoverage
#'
#' Calculates mean allelic coverage (mat+pat) among given replicates.
#'
#' @param df Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)
#' @param reps Optional (default=NA, all replicates), a vector of replicate numbers for which the analysis should be applied
#' @param thr Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis
#' @param thrUP Optional (default=NA), threshold for max gene coverage (default = NA)
#' @param thrType Optional (default = "each", also can be "average" for average coverage on replicates), threshold type
#'
#' @return A table with IDs and calculated mean allelic coverage for given set of replicates
#'
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
