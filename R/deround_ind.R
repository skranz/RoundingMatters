example.deround.ind = function() {
  dat = readRDS("C:/research/comment_phack/MM.Rds")
  dat$num.deci = pmax(num.deci(dat$mu), num.deci(dat$sd))
  dat$s = significand(dat$sd, dat$num.deci)
  dat$sigma = dat$sd
  dat = bind_cols(dat,min.max.z(dat$mu, dat$sigma, dat$num.deci))


  rows = which(dat$s==4 & dat$z==2)
  i = rows[1]
  obs = dat[i,]

  dist = deround.dist.for.obs(dat,obs, min.n=100000, max.rounds=10000)
  unif.dist = sample.uniform.z.deround(nrow(dist),abs(obs$mu), obs$sigma)

  plot.dat = bind_rows(tibble(mode="emp.deround", z=dist$z, alpha=0.5), tibble(mode="uniform.deround", z=unif.dist, alpha=0.5))

  library(ggplot2)
  ggplot(plot.dat, aes(x=z, fill=mode, group=mode)) + geom_histogram(position="identity", alpha=0.45) + theme_bw()

  ggplot(plot.dat, aes(x=z, fill=mode, group=mode)) + geom_density(position="identity", alpha=0.45) + theme_bw()

  ggplot(plot.dat, aes(x=z, fill=mode, group=mode)) +
    geom_density(position="identity", alpha=0.45, kernel="rectangular",bw=0.02) +
    #stat_density(kernel="rectangular",bw=0.01, position="identity",alpha=0.45) +
    theme_bw()


  h = 0.2; z0 = 1.96
  plot.dat %>%
    group_by(mode) %>%
    summarize(
      share.above = mean(z >= 1.96))

}

deround.dist.for.obs = function(dat, obs, s.min=100, min.n=10000, min.rounds = 5, max.rounds=1000, deci.col = "num.deci" ) {
  restore.point("deround.dist.for.obs")

  if (!has.col(dat,"z.min")) {
    dat = bind_cols(dat,min.max.z(dat$mu, dat$sigma, dat[[deci.col]]))
    obs = bind_cols(obs,min.max.z(obs$mu, obs$sigma, obs[[deci.col]]))
  }

  in.range = dat$z.max >= obs$z.min & dat$z.min <= obs$z.max
  d = dat[in.range & dat$s >= s.min,]

  obs.sigma = obs$sigma
  obs.dec = obs[[deci.col]]
  obs.z = abs(obs$z)

  mu = abs(d$mu)
  sigma = d$sigma
  mu.dec = sigma.dec = d[[deci.col]]
  n = nrow(d)

  res.li = vector("list", max.rounds)
  counter = 0
  total.n = 0

  while(counter < max.rounds) {
    counter = counter+1
    # 1. Uniformely de-round mu and sigma
    mu.dr = mu +  runif(n, -0.5,0.5)*10^(-mu.dec)
    sigma.dr = sigma +  runif(n, -0.5,0.5)*10^(-sigma.dec)
    z.dr = mu.dr / sigma.dr

    # 2. Let's draw n uniformely derounded values of obs$sigma
    #    to which we scale sigma.dr

    sigma.scaled = obs.sigma + runif(n, -0.5,0.5)*10^(-obs.dec)

    # 3. Scale mu in the same way
    mu.scaled = mu.dr * sigma.scaled / sigma.dr

    # 4. Round sigma.scaled and mu.scaled with the same number
    #    of digits as obs and compute resulting z
    z.res = round(mu.scaled, obs.dec) / round(sigma.scaled, obs.dec)

    rows = which(z.res == obs.z)

    res.li[[counter]] = cbind(rows, z.dr[rows])
    total.n = total.n + length(rows)
    if (total.n >= min.n & counter >= min.rounds) break
  }

  res = do.call(rbind,res.li)
  colnames(res) = c("ind","z")

  as.data.frame(res)
}


#' Sample derounded z from the uniformely derounded distributon for a given
#' value of mu and sigma
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

