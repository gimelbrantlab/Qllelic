#' MixBetaBinomialFit
#'
#' Fitting beta-binomial mixture distribution of AI.
#'
#' @param initials Initials for EM algm: initials = c(w1, alpha1, alpha2), weight of first component and alphas for both beta-binomial distributions in a mixture
#' @param coverage A number, that reflects higher boundary of the coerage bin
#' @param observations A vector of maternal counts in the bin
#'
#' @return Weight of first component and alphas for both beta-binomial distributions in a mixture, to which algm coincided, plus number of steps.
#' @export
#'
MixBetaBinomialFit <- function(initials, coverage, observations){
  initials_old <- initials
  initials_new <- MixBetaBinomialFitStep(initials_old, coverage, observations)
  n_steps <- 1
  fit_tab <- data.frame(cov = coverage/2, w1 = c(initials_new[1]), a1 = c(initials_new[2]), a2 = c(initials_new[3]))
  while(all(abs(initials_new - initials_old) > 0.001)){
    initials_old <- initials_new
    initials_new <- MixBetaBinomialFitStep(initials_old, coverage, observations)
    if (sum(initials_new < 0) > 0 | sum(is.nan(initials_new)) > 0 | initials_new[3] > 1){
      initials_new <- initials_old
      break
    }
    fit_tab <- rbind(fit_tab, data.frame(cov = coverage/2, w1 = c(initials_new[1]), a1 = c(initials_new[2]), a2 = c(initials_new[3])))
    n_steps = n_steps + 1
  }
  return(list(c(initials_new, n_steps), fit_tab))
}
