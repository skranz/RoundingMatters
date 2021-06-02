#' Density estimates for absolute z-statistics assuming that z-statistics are symmetrically distributed around 0
#'
#' Avoids downward bias at the left hand side where abs(z)=0.
#'
#' @param z vector of z-statistics (or absolute z-statistics)
#' @param at vector of points where density shall be evaluated. If NULL return a function (by calling \code{approxfun}) that allows evaluate the density at arbitrary points.
#' @param bw,adjust,kernel,n,... arguments passed to \code{stats::density}
absz.density <- function(z,at=NULL, bw = 0.1, adjust = 1, kernel = "epanechnikov", n = 1024,...) {
  restore.point("absz.density")
  z = abs(z)
  to = max(z)
  z.vec = c(-abs(z[z>0]),z[z==0],abs(z[z>0]))
  z.dens <- stats::density(z.vec, bw = bw, adjust = adjust,kernel = kernel,from=-to, to=to, n = n,...)

  rows = which(z.dens$x>=0)
  dens.adjust = length(z.dens$x)/length(rows)
  if (is.null(at)) {
    approxfun(x=z.dens$x[rows],y=z.dens$y[rows] *dens.adjust,yleft = z.dens$y[rows[1]]*dens.adjust,yright = 0)
  } else {
    approx(x=z.dens$x[rows],y=z.dens$y[rows] *dens.adjust,xout=at,yleft = z.dens$y[rows[1]]*dens.adjust,yright = 0)$y
  }
}



#' ggplot2 density lines for absolute z-statistics assuming that they are symmetrically distributed around 0
#'
#' Unlike normal [geom_density] or [stat_density] the density estimate does
#' not go artificially decrease at the left bound 0. Note that this function
#' only works nicely if the data starts left with 0. Possibly atoms at z=0
#' should ideally be removed.
#'
#' @param bw The smoothing bandwidth to be used.
#'   If numeric, the standard deviation of the smoothing kernel.
#'   If character, a rule to choose the bandwidth, as listed in
#'   [stats::bw.nrd()].
#' @param adjust A multiplicate bandwidth adjustment. This makes it possible
#'    to adjust the bandwidth while still using the a bandwidth estimator.
#'    For example, `adjust = 1/2` means use half of the default bandwidth.
#' @param kernel Kernel. See list of available kernels in [density()].
#' @param n number of equally spaced points at which the density is to be
#'   estimated, should be a power of two, see [density()] for
#'   details
#' @param trim If `FALSE`, the default, each density is computed on the
#'   full range of the data. If `TRUE`, each density is computed over the
#'   range of that group: this typically means the estimated x values will
#'   not line-up, and hence you won't be able to stack density values.
#'   This parameter only matters if you are displaying multiple densities in
#'   one plot or if you are manually adjusting the scale limits.
#' @section Computed variables:
#' \describe{
#'   \item{density}{density estimate}
#'   \item{count}{density * number of points - useful for stacked density
#'      plots}
#'   \item{scaled}{density estimate, scaled to maximum of 1}
#'   \item{ndensity}{alias for `scaled`, to mirror the syntax of
#'    [`stat_bin()`]}
#' }
#' @export
stat_abszdensity <- function(mapping = NULL, data = NULL,
                         geom = "line", position = "stack",
                         ...,
                         bw = "nrd0",
                         adjust = 1,
                         kernel = "epanechnikov",
                         n = 512,
                         trim = FALSE,
                         na.rm = FALSE,
                         orientation = NA,
                         show.legend = NA,
                         inherit.aes = TRUE) {

  layer(
    data = data,
    mapping = mapping,
    stat = StatAbsZDensity,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      bw = bw,
      adjust = adjust,
      kernel = kernel,
      n = n,
      trim = trim,
      na.rm = na.rm,
      orientation = orientation,
      ...
    )
  )
}

StatAbsZDensity <- ggproto("StatAbsZDensity", Stat,
  required_aes = "x|y",

  default_aes = aes(x = after_stat(density), y = after_stat(density), fill = NA, weight = NULL),

  setup_params = function(data, params) {
    params$flipped_aes <- has_flipped_aes(data, params, main_is_orthogonal = FALSE, main_is_continuous = TRUE)

    has_x <- !(is.null(data$x) && is.null(params$x))
    has_y <- !(is.null(data$y) && is.null(params$y))
    if (!has_x && !has_y) {
      abort("stat_density() requires an x or y aesthetic.")
    }

    params
  },

  extra_params = c("na.rm", "orientation"),

  compute_group = function(data, scales, bw = "nrd0", adjust = 1, kernel = "gaussian",
                           n = 512, trim = FALSE, na.rm = FALSE, flipped_aes = FALSE) {
    restore.point("compute_group_symdensity")

    data <- flip_data(data, flipped_aes)
    if (trim) {
      range <- range(data$x, na.rm = TRUE)
    } else {
      range <- scales[[flipped_names(flipped_aes)$x]]$dimension()
    }

    density <- compute_abszdensity(data$x, data$weight, from = 0,
      to = range[2], bw = bw, adjust = adjust, kernel = kernel, n = n)
    density$flipped_aes <- flipped_aes
    flip_data(density, flipped_aes)
  }

)

compute_abszdensity <- function(x, w, from, to, bw = "nrd0", adjust = 1,
                            kernel = "gaussian", n = 512) {
  nx <- length(x)
  if (is.null(w)) {
    w <- rep(1 / nx, nx)
  } else {
    w <- w / sum(w)
  }
  w = NULL

  # if less than 2 points return data frame of NAs and a warning
  if (nx < 2) {
    warn("Groups with fewer than two data points have been dropped.")
    return(new_data_frame(list(
      x = NA_real_,
      density = NA_real_,
      scaled = NA_real_,
      ndensity = NA_real_,
      count = NA_real_,
      n = NA_integer_
    ), n = 1))
  }

  from = 0
  x = abs(x)
  x.vec = c(-abs(x[x>0]),x[x==0],abs(x[x>0]))
  dens <- stats::density(x.vec, weights = NULL, bw = bw, adjust = adjust, from=-to, to=to,kernel = kernel, n = 2*n)

  rows = which(dens$x>=0)

  ggplot2:::new_data_frame(list(
    x = dens$x[rows],
    density = dens$y[rows] * length(dens$x)/length(rows),
    scaled =  dens$y[rows] / max(dens$y[rows], na.rm = TRUE),
    ndensity = dens$y[rows] / max(dens$y[rows], na.rm = TRUE),
    count =   dens$y[rows] * nx,
    n = nx
  ), n = length(dens$x[rows]))
}
