#' CountsToAI
#'
#' Calculates allelic imbalances from merged counts over given replicates (ai(sum_reps(gene)))
#'
#' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns
#' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df)
#' @param meth An optional parameter for method to use, either sum(m)/sum(p), or sum(m/p) (default = sum(m)/sum(p))
#' @param thr An optional parameter for a threshold on gene coverage (default = NA)
#' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
#' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
#'
#' @return Df with names and mean(mean(m_1,...,m_6))_SNP / mean(mean(m_1+p_1,...,m_6+p6))_SNP
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
