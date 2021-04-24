# Functions for "derounding by simulating rounding" approach


example.dsr = function() {

  dat = readRDS("C:/research/comment_phack/MM_new.Rds")

  h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5)

  dat = dsr.mark.obs(dat,h.seq = h.seq)
  dsr = dsr.ab.df(dat, h.seq=h.seq, z0=1.96)

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

}

#' Create an ab.df for the dsr approach
#'
#' The resulting data frame is required for derounding b simulting rounding (dsr) approach.
#' It contains a row for all considered combinations of z and s and window half-width h in h.seq. The columns share.below and share.above indicate which share of derounded z-statistics are inside the window and either fall below or above the threshold z0, respectively. Note that 1-share.above-share.below is the share of derounded z-statistics that fall outside the considered window.
#' @param dat the data frame that should contain at least the columns z and num.deci (number fo decimal places of mu and sigma, maximum of both)
#' @param h.seq all considered window half-widths
#' @param z0 the significance threshold. Can be a single number or a vector with one element per row of dat.
#' @param min.n how many z values shall be minimally rounded to compute the derounded z-distribution for each observation.
#' @param min.rounds how many repetitions of rounding z-values shall there be at least (even if min.n is already reached).
#' @param verbose Shall some progress information be shown? (This function can take a while).
dsr.ab.df = function(dat,h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5), z0=1.96,min.n=10000,min.rounds=5,  verbose=TRUE) {
  restore.point("dsr.ab.df")
  if (!has.col(dat,"dsr.compute")) {
    stop("Please first call dsr.mark.obs on dat to determine the observations for which the dsr derounding statistics shall be computed.")
  }
  deci.col="num.deci"
  if (!has.col(dat, deci.col)) {
    stop(paste0("Please, put the number of decimal places in the column ", deci.col, ". Usually, you would call\n\ndat$num.deci = pmax(num.deci(dat$mu), num.deci(dat$sd)) "))
  }

  rows = which(dat$dsr.compute)

  if (verbose)
    cat("\nCompute dsr summaries for ", length(rows), " observations.")
  res = bind_rows(lapply(seq_along(rows), function(i) {
    obs = dat[rows[i],]
    if (verbose)
      cat("\n",i," row ", rows[i], " z == ", obs$z, " s == ",obs$s)
    dsr.ab.for.obs(dat, obs, h.seq=h.seq, z0=z0,min.n=min.n, min.rounds=min.rounds, deci.col=deci.col)
  }))
  res
}

dsr.ab.for.obs = function(dat,obs,h.seq=c(0.05,0.075,0.1,0.2,0.3,0.4,0.5), z0=1.96, deci.col="num.deci", ...) {
  restore.point("deround.stat.for.obs")

  dist = dsr.dist.for.obs(dat,obs, deci.col=deci.col, ...)

  h.seq = rev(sort(h.seq))
  nh = length(h.seq)
  num.obs = unique.sources = rep(0, nh)
  share.below = share.above = rep(NA_real_,nh)
  n = nrow(dist)
  i = 0
  for (h in h.seq) {
    i = i+1
    dist = dist[abs(dist$z-z0) <=h,]
    num.obs[i] = nrow(dist)
    unique.sources[i] = n_distinct(dist$ind)
    share.below[i] = sum(dist$z<z0) / n
    share.above[i] = sum(dist$z>=z0) / n
  }
  tibble(mode="dsr",rowid = obs$rowid, z=obs$z, s=obs$s,h=h.seq, z0=z0,  obs=num.obs,dsr.sources=unique.sources, share.below=share.below, share.above=share.above)


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
#' @param no.deround a logical vector indicating columns that shall never be derounded
dsr.mark.obs = function(dat, h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5), z0=1.96, s.max = 100, no.deround = NULL) {
  restore.point("dsr.mark.obs")
  n = nrow(dat)
  dsr.adjust = rep(FALSE,n)

  for (h in h.seq) {
    risk.any = rounding.risks(z = dat$z, s=dat$s, z0=z0, h=h)$risk.any
    dsr.adjust = dsr.adjust | risk.any
  }
  dsr.adjust = dsr.adjust & dat$s <= s.max & dat$sigma > 0
  if (!is.null(no.deround))
    dsr.adjust = dsr.adjust & !no.deround
  dat$dsr.adjust = dsr.adjust
  dupl = duplicated(select(dat, z,s,dsr.adjust))
  dat$dsr.compute = dat$dsr.adjust & !dupl
  dat
}

