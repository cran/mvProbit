\name{summary.mvProbitMargEff}
\alias{summary.mvProbitMargEff}
\alias{print.summary.mvProbitMargEff}

\title{Summarize Marginal Effects of Multivariate Probit Models}

\description{
   These methods prepare and print a statistical summary
   of marginal effects of multivariate probit models.
}

\usage{
\method{summary}{mvProbitMargEff}( object, ... )

\method{print}{summary.mvProbitMargEff}( x, digits = 4, ... )
}

\arguments{
   \item{object}{object of class \code{"mvProbitMargEff"}
      (returned by \code{\link{mvProbitMargEff}} or 
      \code{\link{margEff.mvProbit}}).}
   \item{x}{object of class \code{"summary.mvProbitMargEff"}
      (returned by \code{\link{summary.mvProbitMargEff}}).}
   \item{digits}{positive integer specifiying the minimum number of
      significant digits to be printed
      (passed to \code{\link[stats]{printCoefmat}}).}
   \item{\dots}{currently not used.}
}
\value{
   \code{summary.mvProbitMargEff} returns an object 
   of class \code{"summary.mvProbitMargEff"}
   inheriting from class \code{"matrix"}.
   The returned object is a matrix with four columns,
   where the marginal effects are in the first column,
   their standard errors are in the second column,
   corresponding \sQuote{z values} are in the third column,
   and the resulting \sQuote{P values} are in the last column.
}

\author{Arne Henningsen}

\seealso{\code{\link{mvProbitMargEff}}, 
   \code{\link{margEff.mvProbit}},
   \code{\link{mvProbit}}.
}

\keyword{methods}

