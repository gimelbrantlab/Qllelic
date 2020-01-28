#' A function that tests for allelic imbalance
#'
#' @param x Vector two entries: first is the output of CIs test (T if genes are differentially monoallelically expressed), second is absolute AI value
#'
#' @return character (overall MAE status among clones (by Gendrel classification))
#' @export
#'
#' @examples isMAE_test_CI(c(TRUE, 0.1, 0.5))
isMAE_test_CI <- function(x) {
  thr <- x[3]
  if ((is.na(x[1])) | (is.na(x[2]))) {
    return("nd")
  }
  else {
    if ((x[1]) & ((x[2] >= thr))) {
      return("129_monoallelic")
    }
    else if ((x[1]) & ((x[2] <= (1 - thr)))) {
      return("CAST_monoallelic")
    }
    else if ((x[1]) & ((x[2] >= 0.5))) {
      return("129_biased")
    }
    else if ((x[1]) & ((x[2] < 0.5))) {
      return("CAST_biased")
    }
    else {return("biallelic")}
  }
}
