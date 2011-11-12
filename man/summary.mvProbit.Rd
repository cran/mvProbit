\name{summary.mvProbit}
\alias{summary.mvProbit}
\alias{print.summary.mvProbit}

\title{Summary Results of Multivariate Probit Models}

\description{
   These methods prepare and print summary results
   for multivariate probit models.
}

\usage{
\method{summary}{mvProbit}( object, ... )

\method{print}{summary.mvProbit}( x, digits = 4, ... )
}

\arguments{
   \item{object}{object of class \code{"mvProbit"}
      (returned by \code{\link{mvProbit}}).}
   \item{x}{object of class \code{"summary.mvProbit"}
      (returned by \code{\link{summary.mvProbit}}).}
   \item{digits}{positive integer specifiying the minimum number of
      significant digits to be printed
      (passed to \code{\link[stats]{printCoefmat}}).}
   \item{\dots}{currently not used.}
}
\value{
   \code{summary.mvProbit} returns an object of class \code{"summary.mvProbit"}
   inheriting from class \code{"summary.maxLik"}.
   The returned object contains the same components as objects
   returned by \code{\link[maxLik]{summary.maxLik}} and additionally
   the following components:
   \item{call}{the matched call.}
   \item{start}{the vector of starting values.}
   \item{nDep}{the number of dependent variables.}
   \item{nReg}{the number of explanatory variables (regressors).}
   \item{nObs}{the number of observations.}
}

\author{Arne Henningsen}

\seealso{\code{\link{mvProbit}}}

\keyword{methods}

