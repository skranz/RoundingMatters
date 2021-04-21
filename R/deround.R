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
  res.omit = study.with.derounding(filter(dat, s>=37 | is.na(s)), mode="reported", repl=30)
  res.omit$mode = "omit"

  res.omit.unif = study.with.derounding(filter(dat, s>=37 | is.na(s)), mode="uniform", repl=30)
  res.omit.unif$mode = "omit.unif"


  res.rep = study.with.derounding(dat, mode="reported", repl=30)


  res.uniform = study.with.derounding(dat, mode="uniform", repl=30)
  res.uniform

  z.pdf = make.z.pdf(dat=dat, min.s=100, show.hist=TRUE)
  res.zda = study.with.derounding(dat, mode="zda",z.pdf=z.pdf, repl=30)
  res.zda

  #dat = dsr.mark.obs(dat)
  #ab.sum = dsr.compute.summaries(dat)

  ab.df = readRDS("dsr.Rds")
  res.dsr = study.with.derounding(dat, mode="dsr",ab.df = ab.df, repl=30)

  res = bind_rows(res.rep,res.omit, res.omit.unif, res.uniform, res.zda, res.dsr) %>%
    arrange(h)


  deround.mode = c("reported", "uniform","zda","dsr")
  derounding.repl = 30


  h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5)




  h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5)
}

#' Analysis with derounded z-statistics for different window half-widths around z0
#'
#' This is the main function you will call if you want to perform a publication bias / p-hacking analysis with derounded z-statistics. It allows flexible combinations of how a single derounded z vector is drawn, which statistics are computed for each combination of window h and derounded z-draw and how those statistics are aggregated over multiple replications.
#'
#' @param dat a data frame containing all observations. Each observation is a test from a regression table in some article. It must have  the columns \code{mu} (reported coefficient) and \code{sigma} (reported standard error). The optional column \code{no.deround} can specify rows whose z statistic shall never be derounded. \code{dat} can also have the columns z, num.deci, mu.deci and sigma.deci. If those columns do not exist, they will be computed from mu and sigma.
#' @param h.seq All considered half-window sizes
#' @param window.fun The function that computes for each draw of a derounded z vector and a window h the statistics of interest. Examples are \code{\link{window.t.ci}} (DEFAULT) or \code{\link{window.binom.test}}. Not that our implementation of dsr derounding (or any other derounding using \code{ab.df}) does not draw derounded z, but only creates a logical vector \code{above} indicating which draws are above or below the z0 threshold. This means if you write a custom function, it should essentially work on that vector.
#' @param mode Mode how a single draw of derounded z is computed: "reported", "uniform","zda","dsr" or some custom name (requires ab.df to be defined)
#' @param alt.mode Either "uniform" (DEFAULT) or "reported". Some derounding modes like "zda" and "dsr" cannot be well defined (or are too time-consuming to compute) for observations with many significant digits or outlier z-statistics. \code{alt.mode} specifies how z values shall be selected for those observations.
#' @param z0 The significance threshold for z
#' @param repl Number of replications of each derounding draw.
#' @param aggregate.fun How shall multiple replications be aggregated. Not yet implemented. Currently we always take the medians of each variale returned by window.fun of all replications.
#' @param ab.df Required if \code{mode=="dsr"} or some custom mode. See e.g. \code{\link{dsr.ab.df}}.
#' @param z.pdf Required if \code{mode=="zda"}. Should be generated via \code{\link{make.z.pdf}}.
#' @param max.s Used if \code{mode=="zda"}. Specifies the maximum significand for which zda derounding shall be performed. For observations with larger significand s, uniform derounding will be performed.
#' @param common.deci Shall we assume that mu and sigma are given with the same number of decimal places. If \code{TRUE} (Default) take the column \code{num.deci} i present in \code{dat} or create it as the pairwise maximum of the decimal places of mu and sigma. If \code{FALSE}, either use the columns \code{mu.deci} and \code{sigma.deci} if present in \code{dat} or generate them from mu and sigma.
#' @export
study.with.derounding = function(dat, h.seq=c(0.05,0.075,0.1,0.2,0.3,0.4,0.5),window.fun = window.t.ci, mode=c("reported", "uniform","zda","dsr")[1], alt.mode = c("uniform","reported")[1], make.z.fun=NULL, z0 = ifelse(has.col(dat,"z0"), dat[["z0"]], 1.96), repl=1, aggregate.fun="median",  ab.df = NULL,z.pdf=NULL,max.s = 100, common.deci=TRUE, verbose=TRUE) {
  restore.point("study.with.derounding")

  if (!all(has.col(dat,c("mu","sigma"))))
    stop("Your data set dat needs the columns mu and sigma.")


  if (!alt.mode %in%  c("uniform","reported")) {
    stop("alt.mode must be 'uniform' or 'reported'")
  }
  if(!has.col(dat, "z"))
    dat$z = dat$mu / dat$sigma

  if (common.deci) {
    if (!has.col(dat, "num.deci"))
      dat$num.deci = pmax(num.deci(dat$mu),num.deci(dat$sigma))

    dat$use.mu.deci = dat$use.sigma.deci = dat$num.deci
  } else {
    if (!has.col(dat, "mu.deci")) {
      dat$use.mu.deci = num.deci(dat$mu)
    } else {
      dat$use.mu.deci = dat$mu.deci
    }
    if (!has.col(dat, "sigma.deci")) {
      dat$use.sigma.deci = num.deci(dat$sigma)
    } else {
      dat$use.sigma.deci = dat$sigma.deci
    }
  }

  dat$s = significand(dat$sigma,dat$use.sigma.deci)

  no.deround.rows = integer(0)
  if (has.col(dat,"no.deround"))
    no.deround.rows = which(dat$no.deround)


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
    res = compute.stats.for.all.h(z=dat$z, z0=z0, h.seq=h.seq,window.fun=window.fun, dat=dat)
    res$mode = mode
    res$repl = 1

    return(res)
  }

  if (mode == "zda") {
    if (!has.col(dat,"z.min"))
      dat = bind_cols(dat,min.max.z(dat$mu, dat$sigma, dat$use.mu.deci, dat$use.sigma.deci))

    has.risk = has.rounding.risk(dat=dat,h.seq=h.seq)

    just.uniform = dat$s > max.s | !has.risk

    res = lapply(1:repl, function(r) {
      z = deround.z.density.adjust(z.pdf, dat$mu, dat$sigma, dat$use.mu.deci, dat$use.sigma.deci, just.uniform=just.uniform,verbose = verbose)
      z[no.deround.rows] = dat$z[no.deround.rows]

      compute.stats.for.all.h(z=z, z0=z0, h.seq=h.seq,window.fun=window.fun, dat=dat)
    })
  } else if (mode == "uniform") {
    res = lapply(1:repl, function(r) {
      # TO DO: Do we really need to deround all observations?
      z = deround.z.uniform(dat$mu, dat$sigma, dat$use.mu.deci, dat$use.sigma.deci)
      z[no.deround.rows] = dat$z[no.deround.rows]

      compute.stats.for.all.h(z=z, z0=z0, h.seq=h.seq,window.fun=window.fun, dat=dat)
    })
  } else if (!is.null(make.z.fun)) {
    res = lapply(1:repl, function(r) {
      z = make.z.fun(dat)
      compute.stats.for.all.h(z=z, z0=z0, h.seq=h.seq,window.fun=window.fun, dat=dat)
    })
  } else {
    n = NROW(dat)
    dat.zs = paste0(dat$z,"|",dat$s)
    ab.df$zs = paste0(ab.df$z, "|", ab.df$s)

    h = h.seq[1]

    res = lapply(1:repl, function(r) {
      rand = runif(n,0,1)
      if (alt.mode == "uniform") {
        alt.z = deround.z.uniform(dat$mu, dat$sigma, dat$use.mu.deci, dat$use.sigma.deci)
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
        if (length(no.deround.rows)>0) {
          is.above[no.deround.rows] = dat$z[no.deround.rows] >= z0.vec[no.deround.rows]
          in.window[no.deround.rows] =  abs(dat$z[no.deround.rows]-z0.vec[no.deround.rows]) <= h
        }

        use.rows = which(in.window)
        above = is.above[use.rows]
        inner.res = window.fun(above=above,h=h,dat=dat, rows=use.rows)
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


compute.stats.for.all.h = function(z,z0,h.seq=0.1, window.fun,dat=NULL,...) {
  restore.point("compute.stats.for.all.h")

  z0.h = z0
  res = bind_rows(lapply(h.seq, function(h) {
    rows = which(abs(z-z0) <= h)
    z.h = z[rows]
    if (length(z0)>1) {
      z0.h = z0[rows]
    }
    above = z.h >= z0.h
    window.fun(above=above, z=z.h,z0=z0.h,h=h,dat=dat, rows=rows,...)
  }))
  res
}


#' Draw derounded z assuming missing digits of mu and sigma are uniformly distributed
#'
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
