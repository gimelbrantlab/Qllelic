#' MergeSumCounts
#'
#' Creates a table of sums of maternal and paternal count for given replicates.
#'
#' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
#' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df).
#' @param thr An optional parameter for a threshold on gene coverage (default = NA)
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#'
#' @return Table with names and sums of mat and pat coverages among given replicates.
#' @export
#'
#' @examples MeanCoverage(allelicCountsTable, reps=c(1,2))
#'
MergeSumCounts <- function(df, reps=NA, thr=NA, thrUP=NA, thrType="each"){
  options(stringsAsFactors = FALSE)

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
