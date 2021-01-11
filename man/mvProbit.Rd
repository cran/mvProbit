\name{mvProbit}
\alias{mvProbit}
\alias{print.mvProbit}
\title{Estimation of Multivariate Probit Models}
\description{
   Estimating multivariate probit models by the maximum likelihood method.
   
   WARNING: this function is experimental and extremely
   (perhaps even unusably) slow!
}

\usage{
mvProbit( formula, data, start = NULL, startSigma = NULL, 
   method = "BHHH", finalHessian = "BHHH", 
   algorithm = "GHK", nGHK = 1000, 
   intGrad = TRUE, oneSidedGrad = FALSE, eps = 1e-6, 
   random.seed = 123, ... )

\method{print}{mvProbit}( x, digits = 4, ... )
}

\arguments{
   \item{formula}{a \code{"formula"}:
      a symbolic description of the model
      (currently, all binary outcome variables must have the same regressors).}
   \item{data}{a \code{data.frame} containing the data.}
   \item{start}{an optional numeric vector specifying the starting values 
      for the model coefficients; 
      if argument \code{startSigma} is not specified,
      this vector can also include the correlation coefficients;
      the order of elements is explained in the section \dQuote{details};
      if this argument is not specified,
      coefficients estimated by univariate probit models 
      are used as starting values for the model coefficients.}
   \item{startSigma}{optional starting values for the covariance/correlation matrix 
      of the residuals (must be symmetric and have ones on its diagonal);
      if this argument is not specified
      and the starting values for the correlation coefficients
      are not included in argument \code{start},
      the correlation matrix of the \sQuote{response} residuals, 
      i.e. y - pnorm( X' beta ),
      is used as starting values for sigma.}
   \item{method}{maximisation method / algorithm 
      (see \code{\link[maxLik]{maxLik}}).}
   \item{finalHessian}{Calculation of the final Hessian:
      either \code{FALSE} (no calculation of Hessian), 
      \code{TRUE} (finite-distance calculation of Hessian), or 
      \code{"BHHH"} (calculation based on information equality approach
         and finite-distance gradients, the default).}
   \item{algorithm}{algorithm for computing integrals 
      of the multivariate normal distribution,
      either function \code{GenzBretz()}, \code{Miwa()}, or \code{TVPACK()}
      (see documentation of \code{\link[mvtnorm]{pmvnorm}})
      or character string \code{"GHK"} 
      (see documentation of \code{\link[bayesm]{ghkvec}}).}
   \item{nGHK}{numeric value specifying the number of simulation draws 
      of the GHK algorithm for computing integrals 
      of the multivariate normal distribution.}
   \item{intGrad}{logical. If \code{TRUE}, 
      the computation of the gradients
      with respect to the estimated parameters
      is done internally in function \code{\link{mvProbitLogLik}}
      when it computes the log-likelihood values.
      If the optimization method requires gradients 
      and this argument is \code{FALSE},
      \code{\link[maxLik]{maxLik}} computes the gradients
      by \code{\link[maxLik]{numericGradient}},
      which is usually slower 
      than the calculation in function \code{\link{mvProbitLogLik}}.
      This argument should be set to \code{FALSE}
      if an optimisation algorithm is used that is not based on gradients.}
   \item{oneSidedGrad}{logical. If this argument
      and argument \code{intGrad} are both \code{TRUE}, 
      the gradients of the log-likelihood function
      with respect to the estimated parameters 
      are obtained by one-sided numeric finit-difference differentiation,
      which is faster but less precise 
      than two-sided numeric finit-difference differentiation.}
   \item{eps}{numeric. The step size for the one-sided numeric
      finit-distance differentiation.
      Unfortunately, it is currently not possible to set the step size
      for the two-sided numeric finit-distance differentiation.}
   \item{random.seed}{an integer used to seed R's random number generator;
      this is to ensure replicability 
      when computing (cumulative) probabilities of the multivariate normal distribution
      which is required to calculate the log likelihood values;
      \code{set.seed( random.seed )} is called each time before
      a (cumulative) probability of the multivariate normal distribution 
      is computed;
      defaults to 123.}
   \item{x}{object of class \code{mvProbit} (returned by \code{mvProbit}).}
   \item{digits}{positive integer specifiying the minimum number of
      significant digits to be printed
      (see \code{\link[base]{print.default}}).}
   \item{\dots}{additional arguments to \code{mvProbit} are passed
      to \code{\link[maxLik]{maxLik}} and \code{\link[mvtnorm]{pmvnorm}};
      additional arguments to \code{print.mvProbit} are currently ignored.}
}

\value{
   \code{mvProbit} returns an object of class \code{"mvProbit"}
   inheriting from class \code{"maxLik"}.
   The returned object contains the same components as objects
   returned by \code{\link[maxLik]{maxLik}} and additionally
   the following components:
   \item{call}{the matched call.}
   \item{start}{the vector of starting values.}
   \item{nDep}{the number of dependent variables.}
   \item{nReg}{the number of explanatory variables (regressors).}
   \item{nObs}{the number of observations.}
   \item{dummyVars}{vector of character strings
      indicating the names of explanatory variables
      that contain only zeros and ones or only \code{TRUE} and \code{FALSE}.
      It is \code{NULL}, if no explanatory variable is indentified 
      as a dummy variable.}
}

\details{
   It is possible to specify starting values 
   (a) both for the model coefficients and the correlation coefficients
   (using argument \code{start} alone or arguments \code{start} and \code{startSigma}
   together),
   (b) only for the model coefficients (using argument \code{start} alone), or
   (c) only for the correlation coefficients (using argument \code{startSigma} alone).

   If the model has \eqn{n} dependent variables (equations)
   and \eqn{k} explanatory variables in each equation,
   the order of the starting values in argument \code{start} must be as follows:
   \eqn{b_{1,1}}, \ldots, \eqn{b_{1,k}}, 
   \eqn{b_{2,1}}, \ldots, \eqn{b_{2,k}}, \ldots,
   \eqn{b_{n,1}}, \ldots, \eqn{b_{n,k}},
   where \eqn{b_{i,j}} is the coefficient
   of the \eqn{j}th explanatory variable in the \eqn{i}th equation.
   If argument \code{startSigma} is not specified,
   argument \code{start} can additionally include following elements:
   \eqn{R_{1,2}}, \eqn{R_{1,3}}, \eqn{R_{1,4}}, \ldots, \eqn{R_{1,n}},
   \eqn{R_{2,3}}, \eqn{R_{2,4}}, \ldots, \eqn{R_{2,n}}, \ldots,
   \eqn{R_{n-1,n}},
   where \eqn{R_{i,j}} is the correlation coefficient corresponding to
   the \eqn{i}th and \eqn{j}th equation.

   The \sQuote{state} (or \sQuote{seed}) of R's random number generator 
   is saved at the beginning of the \code{mvProbit} function 
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

\seealso{\code{\link{mvProbitLogLik}},
   \code{\link{mvProbitMargEff}},
   \code{\link[sampleSelection]{probit}},
   \code{\link[stats]{glm}}}

\examples{
## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 50

# generate explanatory variables
xMat <- cbind( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = rnorm( nObs ) )

# model coefficients
beta <- cbind( c(  0.8,  1.2, -0.8 ),
               c( -0.6,  1.0, -1.6 ),
               c(  0.5, -0.6,  1.2 ) )

# covariance matrix of error terms
library( miscTools )
sigma <- symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

# generate dependent variables
yMatLin <- xMat \%*\% beta 
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
colnames( yMat ) <- paste( "y", 1:3, sep = "" )

# estimation (BHHH optimizer and GHK algorithm)
estResult <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   data = as.data.frame( cbind( xMat, yMat ) ), iterlim = 1, nGHK = 50 )
summary( estResult )

# same estimation with user-defined starting values
estResultStart <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), iterlim = 1, nGHK = 50 )
summary( estResultStart )
}

\keyword{models}
\keyword{regression}
