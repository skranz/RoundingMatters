# Functions that compute relevant statistics for a given
# draw of derounded z inside a window around z0


#' theta shall be the probability that a z-statistic is
#' above a threshold z0 in a window with half-width h around z0
#'
#' This function returns for a given vector of z
#' (possibly one derounding draw) statistics:
#' estimate, standard error, 95% CI of theta
#' and number of observations etc
#'
#' This can be used inside Monte Carlo simulations
#'
#' @param z A vector of z-values, possible a derounding draw
#' @param z0 The critical treshold of z that marks an observation significant. This can be a single value, or a vector of the same length as z. This allows to use critical t-values that are adjusted for the sample size of the test an observation refers to
#' @param h The half-width of the considered window around z0. If NULL assume that all observations are already in the correct window.
#' @export
compute.t.stats = function(above=z>=z0,h=NA, ci.level=0.9,z,z0,...) {
  restore.point("compute.t.stats")
  n = length(above)
  theta = mean(above)
  se = sd(above) / sqrt(n)

  # We use t-confidence intervals
  # for large n, t.quant is roughly 1.96 for
  # ci.level = 0.95
  t.quant = abs(qt( (1-ci.level)/2,df = n-1))

  ci.low = theta - t.quant*se
  ci.up = theta + t.quant*se

  as_tibble(list(h=h,theta=theta, se=se, ci.level = ci.level, ci.low=ci.low, ci.up = ci.up, obs=n))

}

compute.binom.test = function(above=z>=z0,h=NA, ci.level=0.95,z,z0,...) {
  restore.point("compute.binom.test")
  n = length(above)
  theta = mean(above)
  se = sd(above) / sqrt(n)


  res = binom.test(x=sum(above), n=length(above), p = 0.5,
           alternative = c("greater"),
           conf.level = 0.95)

  ci.low = res$conf.int[1]
  ci.up = res$conf.int[2]

  as_tibble(list(h=h,theta=theta, p.value = res$p.value, ci.level = ci.level, ci.low=ci.low, ci.up = ci.up, obs=n))

}




