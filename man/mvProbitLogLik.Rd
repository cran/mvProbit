\name{mvProbitLogLik}
\alias{mvProbitLogLik}
\alias{logLik.mvProbit}
\title{Log Likelihood Values for Multivariate Probit Models}
\description{
   Function \code{mvProbitLogLik} calculates log likelihood values 
   of multivariate probit models.

   The \code{logLik} model returns or calculates log likelihood values 
   of multivariate probit models estimated by \code{\link{mvProbit}}.
}

\usage{
mvProbitLogLik( formula, coef, sigma = NULL, data,
   algorithm = "GHK", nGHK = 1000, 
   returnGrad = oneSidedGrad, oneSidedGrad = FALSE, eps = 1e-6, 
   random.seed = 123, ... )

\method{logLik}{mvProbit}( object, coef = NULL, data = NULL, 
   algorithm = NULL, nGHK = NULL, random.seed = NULL, ... )
}

\arguments{
   \item{formula}{a \code{"formula"}:
      a symbolic description of the model
      (currently, all binary outcome variables must have the same regressors).}
   \item{coef}{a numeric vector of the model coefficients;
      if argument \code{sigma} is not specified,
      this vector must also include the correlation coefficients;
      the order of elements is explained in the section \dQuote{details}.}
   \item{sigma}{optional argument for specifying
      the covariance/correlation matrix of the residuals
      (must be symmetric and have ones on its diagonal);
      if this argument is not specified, 
      the correlation coefficients must be specified by argument \code{coef}.}
   \item{data}{a \code{data.frame} containing the data.}
   \item{algorithm}{algorithm for computing integrals 
      of the multivariate normal distribution,
      either function \code{GenzBretz()}, \code{Miwa()}, or \code{TVPACK()}
      (see documentation of \code{\link[mvtnorm]{pmvnorm}})
      or character string \code{"GHK"} 
      (see documentation of \code{\link[bayesm]{ghkvec}}).}
   \item{nGHK}{numeric value specifying the number of simulation draws 
      of the GHK algorithm for computing integrals 
      of the multivariate normal distribution.}
   \item{returnGrad}{logical. If \code{TRUE}, the returned object
      has an attribute \code{"gradient"},
      which is a matrix and provides 
      the gradients of the log-likelihood function
      with respect to all parameters 
      (\code{coef} and upper triangle of \code{sigma})
      evaluated at each observation
      and obtained by (two-sided) numeric finite-difference differentiation.}
   \item{oneSidedGrad}{logical. If \code{TRUE},
      attribute \code{"gradient"} of the returned object
      is obtained by one-sided numeric finite-difference differentiation.}
   \item{eps}{numeric. The step size for the numeric
      finite-distance differentiation.}
   \item{random.seed}{an integer used to seed R's random number generator;
      this is to ensure replicability 
      when computing (cumulative) probabilities of the multivariate normal distribution
      which is required to calculate the log likelihood values;
      \code{set.seed( random.seed )} is called each time before
      a (cumulative) probability of the multivariate normal distribution 
      is computed;
      defaults to 123.}
   \item{object}{an object of class \code{"mvProbit"}
      (returned by \code{\link{mvProbit}}.}
   \item{\dots}{additional arguments are passed
      to \code{\link[mvtnorm]{pmvnorm}}
      when calculating conditional expectations.}
}

\value{
   \code{mvProbitLogLik} returns a vector 
   containing the log likelihood values for each observation.

   If argument \code{returnGrad} is \code{TRUE},
   the vector returned by \code{mvProbitLogLik} 
   has an attribute \code{"gradient"},
   which is a matrix and provides 
   the gradients of the log-likelihood function
   with respect to all parameters 
   (\code{coef} and upper triangle of \code{sigma})
   evaluated at each observation
   and obtained by numeric finite-difference differentiation.

   The \code{logLik} method returns the total log likelihood value
   (sum over all observations) 
   with attribute \code{df} equal to the number of estimated parameters
   (model coefficients and correlation coefficients).
}

\details{
   If the \code{logLik} method is called with \code{object}
   as the only argument,
   it returns the log-likelihood value 
   found in the maximum likelihood estimation.
   If any other argument is not \code{NULL},
   the \code{logLik} method calculates the log-likelihood value
   by calling \code{mvProbitLogLik}.
   All arguments that are \code{NULL},
   are taken from argument \code{object}.

   If the model has \eqn{n} dependent variables (equations)
   and \eqn{k} explanatory variables in each equation,
   the order of the model coefficients in argument \code{coef} must be as follows:
   \eqn{b_{1,1}}, \ldots, \eqn{b_{1,k}}, 
   \eqn{b_{2,1}}, \ldots, \eqn{b_{2,k}}, \ldots,
   \eqn{b_{n,1}}, \ldots, \eqn{b_{n,k}},
   where \eqn{b_{i,j}} is the coefficient
   of the \eqn{j}th explanatory variable in the \eqn{i}th equation.
   If argument \code{sigma} is not specified,
   argument \code{coef} must additionally include following elements:
   \eqn{R_{1,2}}, \eqn{R_{1,3}}, \eqn{R_{1,4}}, \ldots, \eqn{R_{1,n}},
   \eqn{R_{2,3}}, \eqn{R_{2,4}}, \ldots, \eqn{R_{2,n}}, \ldots,
   \eqn{R_{n-1,n}},
   where \eqn{R_{i,j}} is the correlation coefficient corresponding to
   the \eqn{i}th and \eqn{j}th equation.

   The \sQuote{state} (or \sQuote{seed}) of R's random number generator 
   is saved at the beginning of the \code{mvProbitLogLik} function 
   and restored at the end of this function 
   so that this function does \emph{not} affect the generation 
   of random numbers outside this function
   although the random seed is set to argument \code{random.seed}
   and the calculation of the (cumulative) multivariate normal distribution
   uses random numbers.
}

\references{
   Greene, W.H. (1996): 
   \emph{Marginal Effects in the Bivariate Probit Model},
   NYU Working Paper No. EC-96-11. 
   Available at \url{https://www.ssrn.com/abstract=1293106}.
}

\author{Arne Henningsen}

\seealso{\code{\link{mvProbit}},
   \code{\link{mvProbitMargEff}},
   \code{\link[sampleSelection]{probit}},
   \code{\link[stats]{glm}}}

\examples{
## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 10

# generate explanatory variables
xMat <- cbind( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = as.numeric( rnorm( nObs ) > 0 ),
   x3 = rnorm( nObs ),
   x4 = rnorm( nObs ) )

# model coefficients
beta <- cbind( c(  0.8,  1.2, -1.0,  1.4, -0.8 ),
               c( -0.6,  1.0,  0.6, -1.2, -1.6 ),
               c(  0.5, -0.6, -0.7,  1.1,  1.2 ) )

# covariance matrix of error terms
library( miscTools )
sigma <- symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

# generate dependent variables
yMatLin <- xMat \%*\% beta 
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
colnames( yMat ) <- paste( "y", 1:3, sep = "" )

# log likelihood values
myData <- as.data.frame( cbind( xMat, yMat ) )
logLikVal <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = myData )
print( logLikVal )
}

\keyword{models}
\keyword{regression}
