#' MergeSumCounts
#'
#' Creates a table of sums of maternal and paternal alellic counts for given replicates.
#'
#' @param df Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)
#' @param reps Optional (default=NA, all replicates), a vector of replicate numbers for which the analysis should be applied
#' @param thr Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis
#' @param thrUP Optional (default=NA), threshold for max gene coverage (default = NA)
#' @param thrType Optional (default = "each", also can be "average" for average coverage on replicates), threshold type
#'
#' @return A table with IDs and calculated allelic counts for each of the given replicates
#'
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
