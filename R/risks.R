#' Compute thresholds for the significant s of the
#' reported standard deviation such that we can rule-out
#' the errors: misclassification, wrong inclusion, wrong exclusion
#'
#' @param z a vector of z statistics
#' @param z0 significance threshold. Can be a single number or a vector of length z
#' @param h half-width of considered window around z0
#' @returns A data frame with the columns "z", "s.misclass", "s.include", "s.exclude" specifying for each z value the corresponding thresholds.
rounding.risk.s.thresholds = function(z, z0=z0, h=0.2) {
  z.below = z < z0
  in.window = abs(z-z0) <= h

  s.misclass = (1+z0) / (2*abs(z0-z))
  s.include = ifelse( in.window,
    ifelse(z.below,
      (1+z0-h) / (2*abs(z0-h-z)),
      (1+z0+h) / (2*abs(z0+h-z))
    ),
    -Inf
  )
  s.exclude = ifelse( !in.window,
    ifelse(z.below,
      (1+z0-h) / (2*abs(z0-h-z)),
      (1+z0+h) / (2*abs(z0+h-z))
    ),
    -Inf
  )
  data.frame(z,z0,h, z.below, in.window, s.misclass, s.include, s.exclude)
}

#' Assess for observations with reported z-statistic z and
#' a signficand of s for the standard error whether it is at risk of
#' the errors:
#' misclassification, wrong inclusion, wrong exclusion
#'
#' @param z a vector of z statistics
#' @param s vector of corresponding significands of the standard error
#' @param z0 significance threshold. Can be a single number like 1.96 or a vector of length z
#' @param h half-width of considered window around z0
#' @returns A data frame with risk of missclassification information for each observations. We illustrate the columns for the misclassification risk:
#'  "s.misclass" is the threshold for the significand s above which we can rule out misclassification risk
#'  \code{risk.misclass = s < s.misclass} indicates whether the observation is at risk of misclassification
#'   \code{risk.misclass.below = risk.misclass & z < z0} indicates whether the observation is at risk of misclassification and below the significance threshold
#'   the other columns should be self-explainable given this info.
rounding.risks = function(z,s, z0=1.96, h=0.2) {
  thresh.df = rounding.risk.s.thresholds(z,z0,h)

  d = thresh.df %>%
    mutate(
      s = s,
      risk.misclass = s < s.misclass,
      risk.exclude = s < s.exclude,
      risk.include = s < s.include,
      risk.any = risk.misclass | risk.exclude | risk.include,

      risk.misclass.below = risk.misclass & z.below,
      risk.misclass.above = risk.misclass & !z.below,
      risk.exclude.below = risk.exclude & z.below,
      risk.exclude.above = risk.exclude & !z.below,
      risk.include.below = risk.include & z.below,
      risk.include.above = risk.include & !z.below
    )
  d
}

#' Summary statistics for rounding risks for different thresholds
#'
#' @param rr.dat A data frame returned from a call to \code{rounding.risks}.
#' @param s.tresh a vector of considered s thresholds
#' @param long if TRUE (default return results in a long format)
#' @export
rounding.risks.summary = function(rr.dat,s.thresh = 0:100, long = TRUE) {
  restore.point("rounding.risks.summary")
  d = rr.dat

  risks.wide = bind_rows(lapply(0:100, function(s.bar) {
    obs.all = sum(d$in.window, na.rm=TRUE)

    sum.fun = function(x) sum(x, na.rm=TRUE) / obs.all
    res = d %>%
      filter(s>=s.bar, in.window) %>%
      summarise(s.bar=s.bar, obs=sum(in.window),
        across(starts_with("risk."),sum.fun)
      )

    res.exclude = d %>%
      filter(s>=s.bar) %>%
      summarise(across(starts_with("risk.exclude"),sum.fun))

    res[names(res.exclude)] = res.exclude
    res
  }))

  if (!long) return(risks.wide)

  dr = risks.wide[,-c(2:5)] %>%
    tidyr::pivot_longer(cols=-1, names_to = "var", values_to = "risk") %>%
    mutate(
      side = ifelse(endsWith(var,".below"),"below","above"),
      error.type = str.between(var,".","."),
      error.lab = case_when(
        error.type == "misclass" ~ "Misclassified",
        error.type == "include" ~ "Wrongly included",
        error.type == "exclude" ~ "Wrongly excluded",

      ) %>% factor(levels=c("Misclassified","Wrongly included","Wrongly excluded"))
    )
  dr
}



#' Compute minimum and maximum possible values of z given rounded mu and sigma
#'
#' @param mu Vector of reported estimated coefficients
#' @param sigma Vector of reported standard errors
#' @param mu.dec Number of reported decimal digits for mu. By default the maximum of the
#' @param long if TRUE (default return results in a long format)
#' @export
min.max.z = function(mu, sigma, mu.dec=pmax(num.deci(mu),num.deci(sigma)), sigma.dec=mu.dec) {
  mu.min = mu - 0.5*10^(-mu.dec)
  mu.max = mu + 0.5*10^(-mu.dec)

  sigma.min = pmax(sigma - 0.5*10^(-sigma.dec),0)
  sigma.max = sigma + 0.5*10^(-sigma.dec)

  data.frame(
    z.min = mu.min / ifelse(mu.min >= 0, sigma.max, sigma.min),
    z.max = mu.max / ifelse(mu.max >= 0, sigma.min, sigma.max)
  )
}

has.rounding.risk = function(dat=NULL, z=dat$z,s=dat$s, h.seq = c(0.05,0.075,0.1,0.2,0.3,0.4,0.5), z0=1.96) {
  restore.point("dsr.mark.obs")
  n = nrow(dat)
  has.risk = rep(FALSE,n)

  for (h in h.seq) {
    risk.any = rounding.risks(z = dat$z, s=dat$s, z0=z0, h=h)$risk.any
    has.risk = has.risk | risk.any
  }
  has.risk
}

