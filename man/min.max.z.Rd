% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/risks.R
\name{min.max.z}
\alias{min.max.z}
\title{Compute minimum and maximum possible values of z given rounded mu and sigma}
\usage{
\method{min}{max.z}(
  mu,
  sigma,
  mu.dec = pmax(num.deci(mu), num.deci(sigma)),
  sigma.dec = mu.dec
)
}
\arguments{
\item{mu}{Vector of reported estimated coefficients}

\item{sigma}{Vector of reported standard errors}

\item{mu.dec}{Number of reported decimal digits for mu. By default the maximum of the}

\item{long}{if TRUE (default return results in a long format)}
}
\description{
Compute minimum and maximum possible values of z given rounded mu and sigma
}
