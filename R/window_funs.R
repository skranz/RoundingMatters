# Functions that compute relevant statistics for a given
# draw of derounded z inside a window around z0


#' Window function returning estimated probability that a z-statistic is
#' above a threshold z0 in a window with half-width h around z0 and t-test
#' confidence intervals
#'
#' Can be used as argument \code{window.fun} in \code{\link{compute.with.derounding}}
#'
#' @export
window.t.ci = function(above=z>=z0,h=NA, ci.level=0.95,z,z0,...) {
  #restore.point("compute.t.stats")
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

window.binom.test = function(above=z>=z0,h=NA, ci.level=0.95,z,z0,...) {
  #restore.point("compute.binom.test")
  n = length(above)
  theta = mean(above)
  se = sd(above) / sqrt(n)


  res = binom.test(x=sum(above), n=length(above), p = 0.5,
           alternative = c("greater"),
           conf.level = ci.level)

  ci.low = res$conf.int[1]
  ci.up = res$conf.int[2]

  as_tibble(list(h=h,theta=theta, p.value = res$p.value, ci.level = ci.level, ci.low=ci.low, ci.up = ci.up, obs=n))

}




