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
compute.theta.stats = function(z, z0=1.96, h=NULL) {
  if (!is.null(h)) {
    rows = abs(z-z0) <= h
    z = z[rows]
    if (length(z0)>1) {
      z0 = z0[rows]
    }
  }

  n = length(z)

  above.z0 = z > z0
  theta = mean(above.z0)
  se = sd(above.z0) / sqrt(n)

  # We use t-confidence intervals
  # for large n, t.quant is roughly 1.96
  t.quant = abs(qt(0.025,df = n-1))

  ci95.low = theta - t.quant*se
  ci95.up = theta + t.quant*se
  list(theta=theta, se=se, ci.low=ci95.low, ci.up = ci95.up, obs=n, share.z2=sum(z==2)/n)

}

#' Draw derounded z assuming missing digits of mu and sigma are uniformly distributed
#' @param mu.round Reported coefficient, possibly rounded
#' @param sigma.round Reported standard error, possibly rounded.
#' @param mu.dec Number of decimal places mu is reported to. Usually, we would assume that mu and sigma are rounded to the same number of decimal places. Since trailing zeros may not be detected, we set the default \code{mu.dec=pmax(num.deci(mu.round),num.deci(sigma.round))}.
#' @param sigma.dec By default equal to mu.dec.
#' @export
deround.z.uniform = function(mu.round, sigma.round,mu.dec=pmax(num.deci(mu.round),num.deci(sigma.round)), sigma.dec=sigma.mu) {
  mu.deround = mu.round + runif(n, -0.5,0.5)*10^(-mu.dec)
  sigma.deround = sigma.round + runif(n, -0.5,0.5)*10^(-sigma.dec)
  z.deround = mu.deround / sigma.deround
  z.deround
}

#' Draw derounded z assuming missing digits of mu and sigma are uniformly distributed, but adjust for estimated density of z using rejection sampling
#' @param z.pdf An estimated density of the derounded z-statistics (e.g. using only observations with many significant digits) normalized such that its highest values is 1. Best use \code{make.normalized.z.pdf} to create such a normalized pdf from a vector of observed z-statistics.
#' @param mu.round Reported coefficient, possibly rounded
#' @param sigma.round Reported standard error, possibly rounded.
#' @param mu.dec Number of decimal places mu is reported to. Usually, we would assume that mu and sigma are rounded to the same number of decimal places. Since trailing zeros may not be detected, we set the default \code{mu.dec=pmax(num.deci(mu.round),num.deci(sigma.round))}.
#' @param sigma.dec By default equal to mu.dec.
#' @param max.rejection.rounds A limit how often the rejection sampler redraws to avoid an infinite loop.
#' @param verbose If \code{TRUE} cat an r for each resampling draw to see how the function progresses.
#' @export
deround.z.density.adjust = function(z.pdf, mu.round, sigma.round, mu.dec=pmax(num.deci(mu.round),num.deci(sigma.round)), sigma.dec=sigma.mu, max.rejection.rounds = 10000, verbose=TRUE) {
  mu.deround = mu.round + runif(n, -0.5,0.5)*10^(-sigma.dec)
  sigma.deround = sigma.round + runif(n, -0.5,0.5)*10^(-sigma.dec)
  z.deround = mu.deround / sigma.deround

  reject.ind = which(runif(n,0,1) > z.pdf(z.deround))
  counter = 0
  while ( (nr<-length(reject.ind)) > 0) {
    counter = counter+1
    if (counter > max.rejection.rounds) {
      warning(paste0("Stopped rejection draws after specified maximum of ", max.rejection.rounds," rounds."))
      break
    }

    if (verbose) {
      cat("r")
      if (counter %% 50 == 0) cat("\n")
    }

    mu.deround = mu.round[reject.ind] + runif(nr, -0.5,0.5)*10^(-sigma.dec[reject.ind])
    sigma.deround = sigma.round[reject.ind] + runif(nr, -0.5,0.5)*10^(-sigma.dec[reject.ind])
    z.deround[reject.ind] = mu.deround / sigma.deround

    reject.ind = reject.ind[which(runif(nr,0,1) > z.pdf(z.deround[reject.ind]))]
  }
  if (verbose) cat("\n")
  z.deround
}
