% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deround_dsr.R
\name{dsr.mark.obs}
\alias{dsr.mark.obs}
\title{Finds observations in dat for which we shall perform dsr adjustment}
\usage{
dsr.mark.obs(
  dat,
  h.seq = c(0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5),
  z0 = 1.96,
  s.max = 100,
  no.deround = NULL
)
}
\arguments{
\item{dat}{the data frame, must have columns \code{z} and \code{s}}

\item{h.seq}{vector of considered windows half-width. We mark an observation for adjustment if it is at risk of missclassification, wrong inclusion, or wrong exclusion for any considered window size.}

\item{z0}{the signficance threshold for z (default=1.96).}

\item{s.max}{only mark observations for adjustment who have \code{s <= s.max}. Default value is 100.}

\item{no.deround}{a logical vector indicating columns that shall never be derounded}
}
\description{
Adds to \code{dat} the logical columns \code{dsr.adjust} and \code{dsr.compute}. \code{dsr.adjust==TRUE} means that z-statistics of this observation will be adjusted by dsr. The adjustment statistics only depend on the reported z value and significant s of sigma. We thus don't need to compute the distribution for all rows with \code{dsr.adjust==TRUE}. If \code{dsr.compute==TRUE} we shall cmpute the derounded distribution for this observations.
}
