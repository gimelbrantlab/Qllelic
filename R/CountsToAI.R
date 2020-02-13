#' CountsToAI
#'
#' Calculates allelic imbalances from merged counts over given replicates (ai(sum_reps(gene))).
#'
#' @param df Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)
#' @param reps Optional (default=NA, all replicates), a vector of replicate numbers for which the analysis should be applied
#' @param meth Optional (default="mergedToProportion", also can be "meanOfProportions"), method to use, either sum(m)/sum(p) (default), or sum(m/p)
#' @param thr Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis
#' @param thrUP Optional (default=NA), threshold for max gene coverage (default = NA)
#' @param thrType Optional (default = "each", also can be "average" for average coverage on replicates), threshold type
#'
#' @return A table with IDs and calculated AI estimate for given set of replicates
#'
#' @export
#'
#' @examples CountsToAI(allelicCountsTable, reps=c(1,2), thr=10, thrUP=1000)
#'
CountsToAI <- function(df, reps=NA, meth="mergedToProportion", thr=NA, thrUP=NA, thrType="each"){
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

  res_df <- data.frame(df[, 1], AI = p, stringsAsFactors = FALSE)
  names(res_df)[1] <- names(df)[1]

  return(res_df)
}
