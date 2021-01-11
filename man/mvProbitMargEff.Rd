\name{mvProbitMargEff}
\alias{mvProbitExp}
\alias{mvProbitMargEff}
\alias{margEff.mvProbit}
\title{Expectations and Marginal Effects from Multivariate Probit Models}
\description{
   \code{mvProbitExp} calculates expected outcomes 
   from multivariate probit models.

   \code{mvProbitMargEff} calculates marginal effects of the explanatory variables
   on expected outcomes from multivariate probit models.

   The \code{margEff} method for objects of class \code{"mvProbit"}
   is a wrapper function
   that (for the convenience of the user)
   extracts the relevant information from the estimation results
   and then calls \code{mvProbitMargEff}.
}

\usage{
mvProbitExp( formula, coef, sigma = NULL, data,
   cond = FALSE, algorithm = "GHK", nGHK = 1000, random.seed = 123, ... )

mvProbitMargEff( formula, coef, sigma = NULL, vcov = NULL, data,
   cond = FALSE, algorithm = "GHK", nGHK = 1000, 
   eps = 1e-06, dummyVars = NA, addMean = FALSE, 
   returnJacobian = FALSE, random.seed = 123,
   ... )

\method{margEff}{mvProbit}( object, data = eval( object$call$data ),
   cond = FALSE, othDepVar = NULL, dummyVars = object$dummyVars, 
   atMean = FALSE, calcVCov = FALSE, ... )
}

\arguments{
   \item{formula}{a one-sided or two-sided \code{"formula"}:
      a symbolic description of the model
      (currently, all binary outcome variables must have 
      the same explanatory variables).}
   \item{coef}{a numeric vector of the model coefficients;
      if argument \code{sigma} is not specified,
      this vector must also include the correlation coefficients;
      the order of elements is explained in the section \dQuote{details}.}
   \item{sigma}{optional argument for specifying
      the covariance/correlation matrix of the residuals
      (must be symmetric and have ones on its diagonal);
      if this argument is not specified, 
      the correlation coefficients must be specified by argument \code{coef}.}
   \item{vcov}{an optional symmetric matrix 
      specifying the variance-covariance matrix of all coefficients 
      (model coefficients and correlation coefficients);
      if this argument is specified,
      the approximate variance covariance matrices of the marginal effects
      are calculated and returned as an attribute (see below).}
   \item{data}{a \code{data.frame} containing the data.}
   \item{cond}{logical value indicating 
      whether (marginal effects on) conditional expectations (if \code{TRUE})
      or (marginal effects on) unconditional expectations (if \code{FALSE}, default)
      should be returned.}
   \item{algorithm}{algorithm for computing integrals 
      of the multivariate normal distribution,
      either function \code{GenzBretz()}, \code{Miwa()}, or \code{TVPACK()}
      (see documentation of \code{\link[mvtnorm]{pmvnorm}})
      or character string \code{"GHK"} 
      (see documentation of \code{\link[bayesm]{ghkvec}}).}
   \item{nGHK}{numeric value specifying the number of simulation draws 
      of the GHK algorithm for computing integrals 
      of the multivariate normal distribution.}
   \item{eps}{numeric, the difference between the two values 
      of each (numerical) explanatory
      variable that is used for the numerical differentiation.}
   \item{dummyVars}{optional vector containing the names (character strings)
      of explanatory variables
      that should be treated as dummy variables (see section \sQuote{Details}).
      If \code{NA} (the default), dummy variables are detected automatically,
      i.e. all explanatory variables
      which contain only zeros and ones or only \code{TRUE} and \code{FALSE}
      in the data set specified by argument \code{data}
      are treated as dummy variables.
      If \code{NULL},
      no variable is treated as dummy variable.}
   \item{addMean}{logical.
      If \code{TRUE}, the mean of values of all marginal effects are added
      in an additional row at the bottom of the returned data.frame.
      If argument \code{returnJacobian} is \code{TRUE},
      the Jacobian of the mean marginal effects with respect to the coefficients 
      is included in the returned array of the Jacobians
      (in an additional slot at the end of the first dimension).
      If argument \code{vcov} of \code{mvProbitMargEff} is specified 
      or argument \code{calcVCov} of the \code{margEff} method is \code{TRUE},
      the variance covariance matrix of the mean marginal effects
      is included in the returned array of the variance covariance matrices
      (in an additional slot at the end of the first dimension).}
   \item{returnJacobian}{logical.
      If \code{TRUE}, the Jacobian of the marginal effects 
      with respect to the coefficients is returned.}
   \item{random.seed}{an integer used to seed R's random number generator;
      this is to ensure replicability 
      when computing (cumulative) probabilities of the multivariate normal distribution
      which is required to calculate the conditional expectations;
      \code{set.seed( random.seed )} is called each time before
      a (cumulative) probability of the multivariate normal distribution 
      is computed;
      defaults to 123.}
   \item{object}{an object of class \code{"mvProbit"}
      (returned by \code{\link{mvProbit}}.}
   \item{othDepVar}{optional scalar or vector for specifying
      the values of the (other) dependent variables
      when calculating the marginal effects on the conditional expectations.
      If this argument is a scalar (zero or one), 
      it is assumed that all (other) dependent variables have this value
      at all observations.
      If this argument is a vector (filled with zeros or ones)
      with length equal to the number of dependent variables, 
      it is assumed that the vector of dependent variables has these values
      at all observations.
      If this argument is \code{NULL} (the default), 
      the dependent variables are assumed to have the values 
      that these variables have in the data set \code{data}.}
   \item{atMean}{logical. 
      If \code{TRUE}, the marginal effects are calculated
      not at each observation 
      but at the mean values across all observations 
      of the variables in the data set specified by argument \code{data}.}
   \item{calcVCov}{logical. 
      If \code{TRUE}, 
      the approximate variance covariance matrices of the marginal effects
      are calculated and returned as an attribute (see below).}
   \item{\dots}{additional arguments to \code{mvProbitExp} are passed
      to \code{\link[mvtnorm]{pmvnorm}}
      when calculating conditional expectations;
      additional arguments of \code{mvProbitMargEff} are passed
      to \code{mvProbitExp} and possibly further
      to \code{\link[mvtnorm]{pmvnorm}};
      additional arguments of the \code{margEff} method are passed
      to \code{mvProbitMargEff}
      and possibly further to \code{mvProbitMargEff}
      and \code{\link[mvtnorm]{pmvnorm}}.}
}

\value{
   \code{mvProbitExp} returns a data frame 
   containing the expectations of the dependent variables.

   \code{mvProbitMargEff} and the \code{margEff} method return a data frame 
   containing the marginal effects of the explanatory variables
   on the expectations of the dependent variables.

   If argument \code{vcov} of function \code{mvProbitMargEff} is specified
   or argument \code{calcVCov} of the \code{margEff} method is \code{TRUE},
   the returned data frame has an attribute \code{vcov},
   which is a three-dimensional array,
   where the first dimension corresponds to the observation
   and the latter two dimensions span the approximate variance covariance matrix
   of the marginal effects calculated for each observation.

   If argument \code{returnJacobian} of function \code{mvProbitMargEff}
   or method \code{margEff} is set to \code{TRUE},
   the returned data frame has an attribute \code{jacobian},
   which is a three-dimensional array
   that contains the Jacobian matrices of the marginal effects 
   with respect to the coefficients at each observation,
   where the first dimension corresponds to the observations,
   the second dimension corresponds to the marginal effects,
   and the third dimension corresponds to the coefficients.
}

\details{
   When calculating (marginal effects on) unconditional expectations,
   the left-hand side of argument \code{formula} is ignored.
   When calculating (marginal effects on) conditional expectations
   and argument \code{formula} is a one-sided formula
   (i.e. only the right-hand side is specified)
   or argument \code{othDepOne} is \code{TRUE},
   (the marginal effects on) the conditional expectations 
   are calculated based on the assumption
   that all other dependent variables are one.

   The computation of the marginal effects
   of dummy variables 
   (i.e. variables specified in argument \code{dummyVars})
   ignores argument \code{eps} 
   and evaluates the effect of increasing these variables from zero to one.
   The marginal effects of (continuous) variables 
   (i.e. variables not specified in argument \code{dummyVars})
   are calculated by evaluating the effect 
   of increasing these variables from their actual values minus \code{0.5 * eps} 
   to their actual values plus \code{0.5 * eps} (divided by \code{eps}).

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

   If argument \code{vcov} of function \code{mvProbitMargEff} is specified
   or argument \code{calcVCov} of the \code{margEff} method is \code{TRUE},
   the approximate variance covariance matrices of the marginal effects
   are calculated at each observation by using the \sQuote{delta method},
   where the jacobian matrix of the marginal effects 
   with respect to the coefficients is obtained by numerical differentiation.

   The \sQuote{state} (or \sQuote{seed}) of R's random number generator 
   is saved at the beginning of the call to these functions
   and restored at the end
   so that these functions do \emph{not} affect the generation 
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
   \code{\link{mvProbitLogLik}},
   \code{\link[sampleSelection]{probit}},
   \code{\link[stats]{glm}}}

\examples{
## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 10

# generate explanatory variables
xData <- data.frame( 
   const = rep( 1, nObs ),
   x1 = as.numeric( rnorm( nObs ) > 0 ),
   x2 = as.numeric( rnorm( nObs ) > 0 ),
   x3 = rnorm( nObs ),
   x4 = rnorm( nObs ) )

# model coefficients
beta <- c(  0.8,  1.2, -1.0,  1.4, -0.8,
           -0.6,  1.0,  0.6, -1.2, -1.6,
            0.5, -0.6, -0.7,  1.1,  1.2 )

# covariance matrix of error terms
library( miscTools )
sigma <- symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

# unconditional expectations of dependent variables
yExp <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = xData )
print( yExp )

# marginal effects on unconditional expectations of dependent variables
margEffUnc <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = xData )
print( margEffUnc )

# conditional expectations of dependent variables
# (assuming that all other dependent variables are one)
yExpCond <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = beta,
   sigma = sigma, data = xData, cond = TRUE )
print( yExpCond )

# marginal effects on conditional expectations of dependent variables
# (assuming that all other dependent variables are one)
margEffCond <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = beta,
   sigma = sigma, data = xData, cond = TRUE )
print( margEffCond )

# conditional expectations of dependent variables
# (assuming that all other dependent variables are zero)
xData$y1 <- 0
xData$y2 <- 0
xData$y3 <- 0
yExpCond0 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = beta, sigma = sigma, data = xData, cond = TRUE )
print( yExpCond0 )

# marginal effects on conditional expectations of dependent variables
# (assuming that all other dependent variables are zero)
margEffCond0 <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = beta, sigma = sigma, data = xData, cond = TRUE )
print( margEffCond0 )
}

\keyword{models}
\keyword{regression}
