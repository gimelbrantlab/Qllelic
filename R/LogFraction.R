#' LogFraction
#'
#' helper function for computing Gamma
#'
#' @param a_up double()
#' @param a_down double()
#' @param j int()
#'
#' @return log((j + a_up) / (j + a_down))
#'
LogFraction <- function(a_up, a_down, j){
  return(log1p((a_up -1)/(j+1)) - log1p((a_down -1)/(j+1)))
}
