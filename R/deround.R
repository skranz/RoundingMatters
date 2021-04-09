#
#

example = function() {
  setwd("C:/research/comment_phack/")
  dat = readRDS("C:/research/comment_phack/MM.Rds")
  dat$num.deci = pmax(num.deci(dat$mu), num.deci(dat$sd))
  dat$mu = abs(dat$mu)
  dat$sigma = abs(dat$sd)
  dat$s = significand(dat$sigma, dat$num.deci)

  all.dat = dat

  dat = filter(all.dat, method=="DID")

  dat = all.dat
  res.omit = compute.with.derounding(filter(dat, s>=37 | is.na(s)), mode="reported", repl=30)
  res.omit$mode = "omit"

  res.omit.unif = compute.with.derounding(filter(dat, s>=37 | is.na(s)), mode="uniform", repl=30)
  res.omit.unif$mode = "omit.unif"


  res.rep = compute.with.derounding(dat, mode="reported", repl=30)


  res.uniform = compute.with.derounding(dat, mode="uniform", repl=30)
  res.uniform

  z.pdf = make.z.pdf(dat=dat, min.s=100, show.hist=TRUE)
  res.zda = compute.with.derounding(dat, mode="zda",z.pdf=z.pdf, repl=30)
  res.zda

  #dat = dsr.mark.obs(dat)
  #ab.sum = dsr.compute.summaries(dat)

  ab.df = readRDS("dsr.Rds") %>% dsr.to.long() %>% filter(zmode=="dsr")
  res.dsr = compute.with.derounding(dat, mode="dsr",ab.df = ab.df, repl=30)

  res = bind_rows(res.rep,res.omit, res.omit.unif, res.uniform, res.zda, res.dsr) %>%
    arrange(h)


  deround.mode = c("reported", "uniform","zda","dsr")
  derounding.repl = 30


  h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5)




  h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5)
}

compute.with.derounding = function(dat, h.seq=c(0.05,0.075,0.1,0.2,0.3,0.4,0.5), mode=c("reported", "uniform","zda","dsr")[1],fun = compute.binom.test, z0 = ifelse(has.col(dat,"z0"), dat[["z0"]], 1.96), z.pdf=NULL, ab.df = NULL, repl=1, aggregate.fun="median", max.s = 100, alt.mode = c("uniform","reported")[2]) {
  restore.point("compute.with.derounding")

  if (!all(has.col(dat,c("z","mu","sigma"))))
    stop("Your data set dat needs the columns z, mu, sigma.")

  if (!alt.mode %in%  c("uniform","reported")) {
    stop("alt.mode must be 'uniform' or 'reported'")
  }

  if (!has.col(dat, "num.deci"))
    dat$num.deci = pmax(num.deci(dat$mu),num.deci(dat$sigma))

  z0.vec = rep(z0,length.out=NROW(dat))


  if (mode == "zda") {
    if (is.null(z.pdf))
      stop("For mode zda you must provide the argument z.pdf. Generate it with a call to make.z.pdf.")
  } else if (! mode %in% c("reported","uniform")) {
    if (is.null(ab.df)) {
      stop("Unless you have mode 'reported', 'uniform' or 'zda', you must provide the argument ab.df.")
    }
  }

  # For reported we run only once
  if (mode == "reported") {
    res = compute.stats.for.all.h(z=dat$z, z0=z0, h.seq=h.seq,fun=fun, dat=dat)
    res$mode = mode
    res$repl = 1

    return(res)
  }

  if (mode == "zda") {
    if (!has.col(dat,"z.min"))
      dat = bind_cols(dat,min.max.z(dat$mu, dat$sigma, dat$num.deci))

    has.risk = has.rounding.risk(dat=dat,h.seq=h.seq)

    just.uniform = dat$s > max.s | !has.risk

    res = lapply(1:repl, function(r) {
      z = deround.z.density.adjust(z.pdf, dat$mu, dat$sigma, dat$num.deci, just.uniform=just.uniform)
      compute.stats.for.all.h(z=z, z0=z0, h.seq=h.seq,fun=fun, dat=dat)
    })
  } else if (mode == "uniform") {
    res = lapply(1:repl, function(r) {
      # TO DO: Do we really need to deround all observations?
      z = deround.z.uniform(dat$mu, dat$sigma, dat$num.deci)
      compute.stats.for.all.h(z=z, z0=z0, h.seq=h.seq,fun=fun, dat=dat)
    })
  } else {
    n = NROW(dat)
    dat.zs = paste0(dat$z,"|",dat$s)
    ab.df$zs = paste0(ab.df$z, "|", ab.df$s)

    h = h.seq[1]

    res = lapply(1:repl, function(r) {
      rand = runif(n,0,1)
      if (alt.mode == "uniform") {
        alt.z = deround.z.uniform(dat$mu, dat$sigma, dat$num.deci)
      } else if (alt.mode == "reported") {
        alt.z = dat$z
      }
      res.i = bind_rows(lapply(h.seq, function(h) {
        abh.df = ab.df[ab.df$h==h,]
        if (NROW(abh.df)==0) {
          warning(paste0("No ab.df data specified for window half-width h=", h))
        }
        mrows = match(dat.zs, abh.df$zs)

        share.below = abh.df$share.below[mrows]
        share.above = abh.df$share.above[mrows]

        is.above = rand <= share.above
        in.window = rand <= share.below+share.above

        # Observations in dat that are not specified in
        # ab.df
        alt.rows = which(is.na(mrows))
        if (length(alt.rows)>0) {
          in.window[alt.rows] = abs(alt.z[alt.rows]-z0.vec[alt.rows]) <= h
          is.above[alt.rows] = alt.z[alt.rows] >= z0.vec[alt.rows]
        }
        use.rows = which(in.window)
        above = is.above[use.rows]
        inner.res = fun(above=above,h=h,dat=dat, rows=use.rows,...)
        inner.res$obs.alt = sum(in.window[alt.rows], na.rm=TRUE)
        inner.res
      }))
    })
  }

  if (repl >1) {
    if (aggregate.fun=="median") {
      res = bind_rows(res)
      res = res %>%
        group_by(h) %>%
        summarize(
          across(.fns=~median(., na.rm=TRUE))
        )
    }
  } else {
    res = res[[1]]
  }

  res$mode = mode
  res$repl = repl
  res

}


compute.stats.for.all.h = function(z,z0,h.seq=0.1, fun,dat=NULL,...) {
  restore.point("compute.stats.for.all.h")

  z0.h = z0
  res = bind_rows(lapply(h.seq, function(h) {
    rows = which(abs(z-z0) <= h)
    z.h = z[rows]
    if (length(z0)>1) {
      z0.h = z0[rows]
    }
    above = z.h >= z0.h
    fun(above=above, z=z.h,z0=z0.h,h=h,dat=dat, rows=rows,...)
  }))
  res
}


#' Draw derounded z assuming missing digits of mu and sigma are uniformly distributed
#' @param mu Reported coefficient, possibly rounded
#' @param sigma Reported standard error, possibly rounded.
#' @param mu.dec Number of decimal places mu is reported to. Usually, we would assume that mu and sigma are rounded to the same number of decimal places. Since trailing zeros may not be detected, we set the default \code{mu.dec=pmax(num.deci(mu),num.deci(sigma))}.
#' @param sigma.dec By default equal to mu.dec.
#' @export
deround.z.uniform = function(mu, sigma,mu.dec=pmax(num.deci(mu),num.deci(sigma)), sigma.dec=mu.dec) {
  n = length(mu)
  mu.deround = mu + runif(n, -0.5,0.5)*10^(-mu.dec)
  sigma.deround = sigma + runif(n, -0.5,0.5)*10^(-sigma.dec)
  z.deround = mu.deround / sigma.deround
  z.deround
}
