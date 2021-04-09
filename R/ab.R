uniform.ab.df = function(dat,h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5), z0=1.96,n=10000, deci.col="num.deci") {
  restore.point("uniform.ab.df")
  rows = which(dat$dsr.compute)
  res = bind_rows(lapply(seq_along(rows), function(i) {
    obs = dat[rows[i],]
    uniform.ab.for.obs(dat, obs, h.seq=h.seq, z0=z0,n=n, deci.col=deci.col)
  }))
  res
}


uniform.ab.for.obs = function(dat,obs,h.seq=c(0.05,0.075,0.1,0.2,0.3,0.4,0.5), z0=1.96, n = 100000, deci.col="num.deci", ...) {
  restore.point("uniform.ab.for.obs")

  dist = sample.uniform.z.deround(n.unif,abs(obs$mu), obs$sigma, obs[[deci.col]])

  h.seq = rev(sort(h.seq))
  nh = length(h.seq)
  num.obs = rep(0, nh)
  share.below = share.above = rep(NA_real_,nh)
  n = NROW(dist)
  i = 0
  for (h in h.seq) {
    i = i+1
    dist = dist[abs(unif.dist-z0) <= h]
    num.obs[i] = NROW(dist)
    share.below[i] = sum(dist<z0) / n
    share.above[i] = sum(dist>=z0) / n
  }
  tibble(mode="uniform",rowid = obs$rowid, z=obs$z, s=obs$s,h=h.seq, z0=z0,  obs=num.obs,share.below=share.below, share.above=share.above)
}


#' Sample derounded z from the uniformely derounded distributon for a given
#' single value of mu and sigma
#'
#' @param n Number of sample draws
#' @param mu Reported coefficient, possibly rounded
#' @param sigma Reported standard error, possibly rounded.
#' @param mu.dec Number of decimal places mu is reported to. Usually, we would assume that mu and sigma are rounded to the same number of decimal places. Since trailing zeros may not be detected, we set the default \code{mu.dec=pmax(num.deci(mu.round),num.deci(sigma.round))}.
#' @param sigma.dec By default equal to mu.dec.
sample.uniform.z.deround = function(n,mu, sigma, mu.dec=pmax(num.deci(mu),num.deci(sigma)), sigma.dec=mu.dec) {
  mu.deround = mu + runif(n, -0.5,0.5)*10^(-mu.dec)
  sigma.deround = sigma + runif(n, -0.5,0.5)*10^(-sigma.dec)
  z.deround = mu.deround / sigma.deround
  z.deround
}


