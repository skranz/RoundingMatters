# Functions for "derounding by simulating rounding" approach


example.dsr = function() {
  dat = readRDS("C:/research/comment_phack/MM.Rds")
  dat$num.deci = pmax(num.deci(dat$mu), num.deci(dat$sd))
  dat$sigma = abs(dat$sd)
  dat$s = significand(dat$sigma, dat$num.deci)
  dat = bind_cols(dat,min.max.z(dat$mu, dat$sigma, dat$num.deci))

  h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5)

  dat = dsr.mark.obs(dat,h.seq = h.seq)
  sum(dat$dsr.adjust, na.rm=TRUE)
  sum(dat$dsr.compute, na.rm=TRUE)

  dsr = dsr.compute.summaries(dat, h.seq=h.seq)

  saveRDS(dsr,"C:/research/comment_phack/dsr.Rds")

  dsr = readRDS("C:/research/comment_phack/dsr.Rds")

  long = dsr.to.long(dsr)

  long %>%
    group_by(h, zmode) %>%
    summarize(
      share.below = mean(share.below, na.rm=TRUE),
      share.above = mean(share.above, na.rm=TRUE)
    )

  dat = left_join(dat, filter(long,zmode=="dsr") %>% select(z,s,share.below, share.above))


  z0 = 1.96
  sum = dsr %>%
    group_by(h) %>%
    summarize(
      across(.fns=~mean(.,na.rm=TRUE))
    )
  mean(dsr$share.above, na.rm=TRUE)
  mean(dsr$unif.share.above, na.rm=TRUE)

  mean(dsr$share.below, na.rm=TRUE)
  mean(dsr$unif.share.below, na.rm=TRUE)


  dat$dsr.compute = FALSE
  dat$dsr.compute[18] = TRUE
  dsr = dsr.compute.summaries(dat, h.seq=0.5)


  rows = which(dat$s==4 & dat$z==2)
  i = rows[1]
  i = 18
  obs = dat[i,]

  res = dsr.summary.for.obs(dat,obs,h.seq = 0.5, min.n=100000, max.rounds=10000)
  res

  sum = res %>%
    group_by(h) %>%
    summarize(
      share.above,
      share.below,
      excess.above = share.above-unif.share.above+(unif.share.below-share.below)
    )



  dist = deround.dist.for.obs(dat,obs, min.n=100000, max.rounds=10000)
  unif.dist = sample.uniform.z.deround(nrow(dist),abs(obs$mu), obs$sigma)

  plot.dat = bind_rows(tibble(mode="emp.deround", z=dist$z, alpha=0.5), tibble(mode="uniform.deround", z=unif.dist, alpha=0.5))

  library(ggplot2)
  ggplot(plot.dat, aes(x=z, fill=mode, group=mode)) +
    geom_density(position="identity", alpha=0.45) +
    theme_bw()


  h = 0.1; z0 = 1.96
  plot.dat %>%
    filter(abs(z-z0)<=h) %>%
    group_by(mode) %>%
    summarize(
      share.above = mean(z >= z0)
    )




}

dsr.to.long = function(dsr) {
  dsr_base = select(dsr, rowid, z, s,h, z0)

  bind_rows(
    bind_cols(dsr_base, transmute(dsr, zmode="dsr",  share.below=dsr.share.below, share.above=dsr.share.above, num.z=dsr.obs, num.sources=dsr.sources,)),
    bind_cols(dsr_base, transmute(dsr, zmode="unif", share.below=unif.share.below, share.above=unif.share.above, num.z=NA, num.sources=NA)),
    bind_cols(dsr_base, transmute(dsr, zmode="reported",  share.below=reported.share.below, share.above=reported.share.above,num.z=NA, num.sources=NA))
  ) %>%
  arrange(rowid,h)

}

dsr.compute.summaries = function(dat,h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5), z0=1.96,n.unif = 100000,min.n=10000,min.rounds=5, deci.col="num.deci", verbose=TRUE, long=TRUE) {
  restore.point("dsr.compute.summaries")
  if (!has.col(dat,"dsr.compute")) {
    stop("Please first call dsr.mark.obs on dat to determine the observations for which the dsr derounding statistics shall be computed.")
  }
  if (!has.col(dat, deci.col)) {
    stop(paste0("Either put the number of decimal digits in the column ", deci.col, " or specify another column in the argument deci.col."))
  }

  rows = which(dat$dsr.compute)

  if (verbose)
    cat("\nCompute dsr summaries for ", length(rows), " observations.")
  res = bind_rows(lapply(seq_along(rows), function(i) {
    obs = dat[rows[i],]
    if (verbose)
      cat("\n",i," row ", rows[i], " z == ", obs$z, " s == ",obs$s)
    dsr.summary.for.obs(dat, obs, h.seq=h.seq, z0=z0, n.unif=n.unif,min.n=min.n, min.rounds=min.rounds, deci.col=deci.col)
  }))

  res = res %>%
    mutate(
      reported.share.below = z < z0 & abs(z-z0)<=h,
      reported.share.above = z >= z0 & abs(z-z0)<=h
    )

  if (long) return(dsr.to.long(res))
  res
}

dsr.summary.for.obs = function(dat,obs,h.seq=c(0.05,0.075,0.1,0.2,0.3,0.4,0.5), z0=1.96, n.unif = 100000, deci.col="num.deci", ...) {
  restore.point("deround.stat.for.obs")

  dist = dsr.dist.for.obs(dat,obs, deci.col=deci.col, ...)
  unif.dist = sample.uniform.z.deround(n.unif,abs(obs$mu), obs$sigma, obs[[deci.col]])

  h.seq = rev(sort(h.seq))
  nh = length(h.seq)
  num.obs = unique.sources = rep(0, nh)
  share.below = share.above = unif.share.below = unif.share.above = rep(NA_real_,nh)
  n = nrow(dist)
  i = 0
  for (h in h.seq) {
    i = i+1
    dist = dist[abs(dist$z-z0) <=h,]
    unif.dist = unif.dist[abs(unif.dist-z0) <= h]
    num.obs[i] = nrow(dist)
    unique.sources[i] = n_distinct(dist$ind)
    share.below[i] = sum(dist$z<z0) / n
    share.above[i] = sum(dist$z>=z0) / n

    unif.share.below[i] = sum(unif.dist<z0) / n.unif
    unif.share.above[i] = sum(unif.dist>=z0) / n.unif
  }
  tibble(rowid = obs$rowid, z=obs$z, s=obs$s,h=h.seq, z0=z0, dsr.obs=num.obs,dsr.sources=unique.sources, dsr.share.below=share.below, dsr.share.above=share.above, unif.share.below, unif.share.above)


}


dsr.dist.for.obs = function(dat, obs, s.min=100, min.n=10000, min.rounds = 5, max.rounds=5000, deci.col = "num.deci" ) {
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


#' Finds observations in dat for which we shall perform dsr adjustment
#'
#' Adds to \code{dat} the logical columns \code{dsr.adjust} and \code{dsr.compute}. \code{dsr.adjust==TRUE} means that z-statistics of this observation will be adjusted by dsr. The adjustment statistics only depend on the reported z value and significant s of sigma. We thus don't need to compute the distribution for all rows with \code{dsr.adjust==TRUE}. If \code{dsr.compute==TRUE} we shall cmpute the derounded distribution for this observations.
#' @param dat the data frame, must have columns \code{z} and \code{s}
#' @param h.seq vector of considered windows half-width. We mark an observation for adjustment if it is at risk of missclassification, wrong inclusion, or wrong exclusion for any considered window size.
#' @param z0 the signficance threshold for z (default=1.96).
#' @param s.max only mark observations for adjustment who have \code{s <= s.max}. Default value is 100.
dsr.mark.obs = function(dat, h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5), z0=1.96, s.max = 100) {
  restore.point("dsr.mark.obs")
  n = nrow(dat)
  dsr.adjust = rep(FALSE,n)

  for (h in h.seq) {
    risk.any = rounding.risks(z = dat$z, s=dat$s, z0=z0, h=h)$risk.any
    dsr.adjust = dsr.adjust | risk.any
  }
  dsr.adjust = dsr.adjust & dat$s <= s.max & dat$sigma > 0
  dat$dsr.adjust = dsr.adjust
  dupl = duplicated(select(dat, z,s,dsr.adjust))
  dat$dsr.compute = dat$dsr.adjust & !dupl
  dat
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

