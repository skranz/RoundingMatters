#' Perform kernel estimates of two densities of absolute z-statistics
#' and their ratio.
#'
#' Add by default bootstrap standard errors and confidence intervals
#'
#' @param z.num observed z-statistics forming numerator density
#' @param z.denom observed z-statistics forming denominator density
#' @param at position where density shall be evaluated
#' @param bootstrap if TRUE add bootstrap SE and CI for all measures
#' @param B number of bootstrap repetitions
#' @param ci.level Confidence level. Default 0.95.
#' @param weights.num weights for z.num (optional)
#' @param weights.denom weights for z.denom (optional)
#' @param ... arguments for absz.density
absz.density.ratio = function(z.num, z.denom, at, bootstrap=TRUE, B=1000,ci.level=.95, bw = 0.1, kernel = "epanechnikov",  return.as = c("long","wide")[1],weights.num=NULL, weights.denom=NULL, ...) {

  z.num = abs(z.num)
  z.denom = abs(z.denom)
  density.num = absz.density(z.num, at=at, bw=bw, kernel = kernel, weights=weights.num, ...)
  density.denom = absz.density(z.denom, at=at, bw=bw, kernel = kernel, weights=weights.denom, ...)
  ratio = density.num / density.denom
  if (return.as == "wide" & !bootstrap) {
    res = quick.df(at = at, density.num=density.num, density.denom = density.denom, ratio = ratio)
  } else {
    N = length(at)
    res = bind_rows(
      quick.df(at=at, var=rep("density.num",N), estimate=density.num),
      quick.df(at=at, var=rep("density.denom",N), estimate=density.denom),
      quick.df(at=at, var=rep("ratio",N), estimate=ratio),
      quick.df(at=at, var=rep("logratio",N), estimate=log(ratio))
    )
  }
  if (!bootstrap) return(res)

  bdf = bind_rows(lapply(1:B, function(i) {
    zb.num = sample(z.num, NROW(z.num),replace = TRUE)
    zb.denom = sample(z.denom, NROW(z.denom),replace = TRUE)
    res = absz.density.ratio(zb.num, zb.denom,at = at, bootstrap=FALSE,bw=bw, kernel=kernel, return.as="long",... )
    res
  }))
  #restore.point("absz.density.ratio")

  sdf = bdf %>%
    group_by(var, at) %>%
    summarize(
      se = sd(estimate),
      ci.lower = quantile(estimate, (1-ci.level)/2),
      ci.upper = quantile(estimate, 1-(1-ci.level)/2)
    )

  left_join(res, sdf, by=c("var","at"))


}
