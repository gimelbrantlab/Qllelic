#' NameColumns
#'
#' Helper function to quickly rename columns in geneCountTab dataframe
#'
#' @param exp_n Experiment number
#' @param rep_n Number of replicates for the experiment
#'
#' @return Vector with names
#' @export
#'
#' @examples colnames(allelicCountsTable)[2:9] <- c(NameColumns(1,2), NameColumns(2,2))
#'
NameColumns <- function(exp_n, rep_n)  {
  paste0("exp", rep(exp_n, 2*rep_n), "_rep", rep(1:rep_n, each = 2), "_", rep(c("ref", "alt"), rep_n))
}
