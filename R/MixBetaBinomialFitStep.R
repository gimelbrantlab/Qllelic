#' MixBetaBinomialFitStep
#'
#' One step of fitting beta-binomial mixture distribution of AI in particular bin.
#'
#' @param initials_old Initials for EM step: initials = c(w1, alpha1, alpha2), weight of first component and alphas for both beta-binomial distributions in a mixture
#' @param coverage A number, that represents the coverage bin
#' @param observations A vector of "maternal counts" in the bin
#'
#' @return Re-fitted initials for next EM step.
#' @export
#'
MixBetaBinomialFitStep <- function(initials_old, coverage, observations){
  w <- c(initials_old[1], 1-initials_old[1])
  alpha <- initials_old[2:3]

  # EM step 1:
  gamma_n_k <- sapply(1:2, function(k){
    a_k <- alpha[k]
    a_notk <- alpha[2-(k+1)%%2]
    w_k <- w[k]
    w_notk <- w[2-(k+1)%%2]

    gamma_k <- sapply(observations, function(xn){
      if(xn == 0){
        logFractionsProd <- sum(sapply(0:(coverage-xn-1), function(j){LogFraction(a_notk, a_k, j)})) +
          sum(sapply(0:(coverage-1), function(j){LogFraction(2*a_k, 2*a_notk, j)}))
      } else if(xn == coverage) {
        logFractionsProd <- sum(sapply(0:(xn-1), function(j){LogFraction(a_notk, a_k, j)})) +
          sum(sapply(0:(coverage-1), function(j){LogFraction(2*a_k, 2*a_notk, j)}))
      } else {
        logFractionsProd <- sum(sapply(0:(xn-1), function(j){LogFraction(a_notk, a_k, j)})) +
          sum(sapply(0:(coverage-xn-1), function(j){LogFraction(a_notk, a_k, j)})) +
          sum(sapply(0:(coverage-1), function(j){LogFraction(2*a_k, 2*a_notk, j)}))
      }
      result <- 1 / (1 + w_notk/w_k * exp(logFractionsProd))
      return(result)
    })
    return(gamma_k)
  })

  # EM step 2:
  sum_gamma_k <- colSums(gamma_n_k)

  w_new <- sum_gamma_k / length(observations)
  mean_new <- c(0.5*coverage, 0.5*coverage)
  var_new <- sapply(1:2, function(k){
    sum(gamma_n_k[, k] * (observations - mean_new[k])**2) / sum_gamma_k[k]
  })
  alpha_new <- sapply(1:2, function(k){
    (coverage**2 - 4 * var_new[k]) / (8 * var_new[k] - 2 * coverage)
  })

  initials_new <- c(w_new[1], alpha_new)

  return(as.numeric(initials_new))
}
