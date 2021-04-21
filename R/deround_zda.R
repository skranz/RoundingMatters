# Derounding using z-density adjustment


#' Compute a normalized pdf from a vector of z-statistics
#'
#' The PDF is normalized such that the point of highest density is 1
#'
#' @param z a vector of z statistics. Usually, you would select
#'          all values from dat whose mu and sigma have sufficiently
#'          many significant digits
#' @param ... other parameters passed to \code{\link{stats::density}}
make.z.pdf = function(z, bw=0.05, kernel="gaussian", n=512,dat, min.s=100,z.min=0, z.max=5,show.hist=FALSE,...) {
  if (missing(z) & !missing(dat)) {
    rows = which(dat$s >= min.s)
    z = dat$z[rows]
  }
  rows = which(z >= z.min & z <= z.max)
  z = z[rows]
  if (show.hist) hist(z)
  z.dens = density(z,bw=bw, kernel=kernel, n=n,...)
  approxfun(x=z.dens$x,y=z.dens$y / max(z.dens$y),yleft = 0,yright = 0)
}


#' Draw derounded z assuming missing digits of mu and sigma are uniformly distributed, but adjust for estimated density of z using rejection sampling
#' @param z.pdf An estimated density of the derounded z-statistics (e.g. using only observations with many significant digits) normalized such that its highest values is 1. Best use \code{make.z.pdf} to create such a normalized pdf from a vector of observed z-statistics.
#' @param mu Reported coefficient, possibly rounded
#' @param sigma Reported standard error, possibly rounded.
#' @param mu.dec Number of decimal places mu is reported to. Usually, we would assume that mu and sigma are rounded to the same number of decimal places. Since trailing zeros may not be detected, we set the default \code{mu.dec=pmax(num.deci(mu),num.deci(sigma))}.
#' @param sigma.dec By default equal to mu.dec.
#' @param max.rejection.rounds A limit how often the rejection sampler redraws to avoid an infinite loop.
#' @param verbose If \code{TRUE} cat an r for each resampling draw to see how the function progresses.
#' @export
deround.z.density.adjust = function(z.pdf, mu, sigma, mu.dec=pmax(num.deci(mu),num.deci(sigma)), sigma.dec=mu.dec, max.rejection.rounds = 10000, verbose=TRUE, just.uniform=rep(FALSE, length(mu)), z.min=0, z.max=5) {
  restore.point("sjsjlksjdkjlksj")
  n = length(mu)
  mu = abs(mu); sigma = abs(sigma)
  z = mu / sigma
  mu.deround = mu + runif(n, -0.5,0.5)*10^(-mu.dec)
  sigma.deround = sigma + runif(n, -0.5,0.5)*10^(-sigma.dec)
  z.deround = mu.deround / sigma.deround

  reject.ind = which( (runif(n,0,1) > z.pdf(z.deround)) & !just.uniform & z >= z.min & z <= z.max)
  counter = 0
  while ( (nr<-length(reject.ind)) > 0) {
    counter = counter+1
    if (counter > max.rejection.rounds) {
      warning(paste0("Stopped rejection draws after specified maximum of ", max.rejection.rounds," rounds."))
      break
    }

    if (verbose) {
      cat(nr,"")
      if (counter %% 20 == 0) cat("\n")
    }

    mu.deround = mu[reject.ind] + runif(nr, -0.5,0.5)*10^(-mu.dec[reject.ind])
    sigma.deround = sigma[reject.ind] + runif(nr, -0.5,0.5)*10^(-sigma.dec[reject.ind])
    z.deround[reject.ind] = mu.deround / sigma.deround

    reject.ind = reject.ind[which(runif(nr,0,1) > z.pdf(z.deround[reject.ind]))]
  }
  if (verbose) cat("\n")
  z.deround
}
