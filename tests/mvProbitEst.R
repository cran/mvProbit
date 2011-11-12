library( "mvProbit" )

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
sigma <- miscTools::symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

# generate dependent variables
yMatLin <- xMat %*% beta 
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
colnames( yMat ) <- paste( "y", 1:3, sep = "" )

# estimation with the BHHH algorithm, two-sided gradients
estResultBHHH <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma,
   data = as.data.frame( cbind( xMat, yMat ) ), tol = 0.5,
   algorithm = GenzBretz() )
print( estResultBHHH )
summary( estResultBHHH )
logLik( estResultBHHH )
estResultBHHHA <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta, sigma[ lower.tri( sigma ) ] ),
   data = as.data.frame( cbind( xMat, yMat ) ), tol = 0.5,
   algorithm = GenzBretz() )
all.equal( estResultBHHH, estResultBHHHA )

# estimation with the BHHH algorithm, one-sided gradients
estResultBHHH1 <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma,
   data = as.data.frame( cbind( xMat, yMat ) ), tol = 0.5,
   algorithm = GenzBretz(), oneSidedGrad = TRUE )
print( estResultBHHH1 )
summary( estResultBHHH1 )
logLik( estResultBHHH1 )
all.equal( estResultBHHH, estResultBHHH1 )

# estimation with the BFGS algorithm, two-sided gradients
estResultBFGS <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, method = "BFGS",
   data = as.data.frame( cbind( xMat, yMat ) ), 
   tol = 0.5, algorithm = GenzBretz() )
print( estResultBFGS )
summary( estResultBFGS )
logLik( estResultBFGS )

# estimation with the BFGS algorithm, one-sided gradients
estResultBFGS1 <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, method = "BFGS", 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   tol = 0.5, algorithm = GenzBretz(), oneSidedGrad = TRUE )
print( estResultBFGS1 )
summary( estResultBFGS1 )
logLik( estResultBFGS1 )
all.equal( estResultBFGS, estResultBFGS1 )

# estimation with the BFGS algorithm, one-sided gradients, no starting values
estResultBFGS1a <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   data = as.data.frame( cbind( xMat, yMat ) ), method = "BFGS",
   tol = 0.5, algorithm = GenzBretz(), oneSidedGrad = TRUE )
print( estResultBFGS1a )
summary( estResultBFGS1a )
logLik( estResultBFGS1a )

# estimation with the BFGS algorithm, one-sided gradients, no starting values for beta
estResultBFGS1b <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   startSigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "BFGS", tol = 0.5, algorithm = GenzBretz(), oneSidedGrad = TRUE )
print( estResultBFGS1b )
summary( estResultBFGS1b )
logLik( estResultBFGS1b )

# estimation with the BFGS algorithm, one-sided gradients, no starting values for sigma
estResultBFGS1s <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "BFGS", tol = 0.5, algorithm = GenzBretz(), oneSidedGrad = TRUE )
print( estResultBFGS1s )
summary( estResultBFGS1s )
logLik( estResultBFGS1s )

# estimation with the BFGS algorithm, Miwa algorithm for obtaining integrals
estResultBFGSm <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), method = "BFGS",
   tol = 0.5, algorithm = Miwa( steps = 64 ) )
print( estResultBFGSm )
summary( estResultBFGSm )
logLik( estResultBFGSm )
all.equal( estResultBFGS, estResultBFGSm )

# estimation with the BFGS algorithm, GHK algorithm for obtaining integrals
estResultBFGSg <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), method = "BFGS",
   tol = 0.5 )
print( estResultBFGSg )
summary( estResultBFGSg )
logLik( estResultBFGSg )
all.equal( estResultBFGS, estResultBFGSg )
all.equal( estResultBFGSm, estResultBFGSg )

# estimation with the Nelder-Mead algorithm
estResultNM <- mvProbit( cbind( y1, y2, y3 ) ~ x1 + x2,
   start = c( beta ), startSigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) ), 
   method = "NM", reltol = 0.05, algorithm = GenzBretz() )
print( estResultNM )
summary( estResultNM )
logLik( estResultNM )


## testing the logLik method
# argument 'coef'
all.equal( logLik( estResultBHHH ), 
   logLik( estResultBHHH, coef = coef( estResultBHHH ) ) )
logLik( estResultBHHH ) -
   logLik( estResultBHHH, coef = coef( estResultBHHH ) * 0.999 )

# argument 'data'
all.equal( logLik( estResultBHHH ), 
   logLik( estResultBHHH, data = as.data.frame( cbind( xMat, yMat ) ) ) )
logLik( estResultBHHH ) -
   logLik( estResultBHHH, data = as.data.frame( cbind( xMat * 0.999, yMat ) ) )

# argument 'algorithm'
all.equal( logLik( estResultBHHH ), 
   logLik( estResultBHHH, algorithm = GenzBretz() ) )
logLik( estResultBHHH ) -
   logLik( estResultBHHH, algorithm = Miwa() )

# argument 'nGHK'
all.equal( logLik( estResultBHHH ), 
   logLik( estResultBHHH, nGHK = 5555 ) )
logLik( estResultBHHH, algorithm = "GHK" ) -
   logLik( estResultBHHH, algorithm = "GHK", nGHK = 2000 )

# argument 'random.seed'
all.equal( logLik( estResultBHHH ), 
   logLik( estResultBHHH, random.seed = 123 ) )
logLik( estResultBHHH ) -
   logLik( estResultBHHH, random.seed = 1234 )


# marginal effects based on estimated coefficients with covariance matrix
# unconditional marginal effects (with Jacobian)
margEffUnc <- margEff( estResultBFGS, calcVCov = TRUE, returnJacobian = TRUE )
print( margEffUnc )
print( attr( margEffUnc, "vcov" )[ 1:5, , ] )
print( drop( attr( margEffUnc, "vcov" )[ nObs, , ] ) )
print( attr( margEffUnc, "jacobian" )[ 1:5, , ] )
print( drop( attr( margEffUnc, "jacobian" )[ nObs, , ] ) )
summary( margEffUnc )
# now with explicitly specifying dummy variables
margEffUncD <- margEff( estResultBFGS, calcVCov = TRUE, returnJacobian = TRUE,
   dummyVars = c( "x1" ) )
all.equal( margEffUncD, margEffUnc )
# now with seemingly no dummy variables
margEffUncD0 <- margEff( estResultBFGS, calcVCov = TRUE, returnJacobian = TRUE,
   dummyVars = NULL )
summary( margEffUncD0 )
# now with seemingly only dummy variables
margEffUncDA <- margEff( estResultBFGS, calcVCov = TRUE, returnJacobian = TRUE,
   dummyVars = c( "x1", "x2" ) )
summary( margEffUncDA )
# now with returned Jacobian but without variance covariance matrix
margEffUncJac <- margEff( estResultBFGS, returnJacobian = TRUE )
all.equal( attr( margEffUncJac, "jacobian" ), attr( margEffUnc, "jacobian" ) )
# now including mean values of the marginal effects
margEffUncM <- margEff( estResultBFGS, calcVCov = TRUE, returnJacobian = TRUE,
   addMean = TRUE )
all.equal( margEffUnc, margEffUncM[ -(nObs+1), ], check.attributes = FALSE )
print( margEffUncM[ nObs:(nObs+1), ] )
all.equal( attr( margEffUnc, "vcov" ), attr( margEffUncM, "vcov" )[ 1:nObs, , ] )
print( attr( margEffUncM, "vcov" )[ nObs:(nObs+1), , ] )
print( drop( attr( margEffUncM, "vcov" )[ nObs+1, , ] ) )
all.equal( attr( margEffUnc, "jacobian" ), attr( margEffUncM, "jacobian" )[ 1:nObs, , ] )
print( attr( margEffUncM, "jacobian" )[ nObs:(nObs+1), , ] )
print( drop( attr( margEffUncM, "jacobian" )[ nObs+1, , ] ) )
all.equal( summary( margEffUnc )[ , ], 
   summary( margEffUncM )[ 1:( 6 * nObs ), ], check.attributes = FALSE )
printCoefmat( summary( margEffUncM )[ -( 1:( 6 * ( nObs - 1 ) ) ), ] )
# now at mean values of explanatory variables
margEffUncMean <- margEff( estResultBFGS, calcVCov = TRUE, 
   data = as.data.frame( t( colMeans( xMat ) ) ) )
summary( margEffUncMean )
# now with argument 'atMean'
margEffUncMeanA <- margEff( estResultBFGS, calcVCov = TRUE, atMean = TRUE )
all.equal( margEffUncMeanA, margEffUncMean )

# conditional marginal effects
# (assuming that all other dependent variables are as observed)
margEffCondObs <- margEff( estResultBFGS, cond = TRUE, algorithm = GenzBretz() )
print( margEffCondObs )

# conditional marginal effects with covariance matrix at sample mean
# (assuming that all other dependent variables are at there modal values)
# (with Jacobian)
margEffCondObsCov <- margEff( estResultBFGS, cond = TRUE,
   atMean = TRUE, othDepVar = c( colMedians( yMat * 1 ) ), 
   algorithm = GenzBretz(), calcVCov = TRUE, returnJacobian = TRUE )
print( margEffCondObsCov )
print( attr( margEffCondObsCov, "vcov" ) )
print( drop( attr( margEffCondObsCov, "vcov" ) ) )
print( attr( margEffCondObsCov, "jacobian" ) )
print( drop( attr( margEffCondObsCov, "jacobian" ) ) )
summary( margEffCondObsCov )
# now with explicitly specifying dummy variables
margEffCondObsCovD <- margEff( estResultBFGS, cond = TRUE,
   atMean = TRUE, othDepVar = c( colMedians( yMat * 1 ) ),
   algorithm = GenzBretz(),
   calcVCov = TRUE, returnJacobian = TRUE, dummyVars = c( "x1" ) )
all.equal( margEffCondObsCovD, margEffCondObsCov )
# now with seemingly no dummy variables
margEffCondObsCovD0 <- margEff( estResultBFGS, cond = TRUE,
   atMean = TRUE, othDepVar = c( colMedians( yMat * 1 ) ), 
   algorithm = GenzBretz(),
   calcVCov = TRUE, returnJacobian = TRUE, dummyVars = NULL )
summary( margEffCondObsCovD0 )
# now with seemingly only dummy variables
margEffCondObsCovDA <- margEff( estResultBFGS, cond = TRUE,
   atMean = TRUE, othDepVar = c( colMedians( yMat * 1 ) ), 
   algorithm = GenzBretz(),
   calcVCov = TRUE, returnJacobian = TRUE, dummyVars = c( "x1", "x2" ) )
summary( margEffCondObsCovDA )

# conditional marginal effects
# (assuming that all other dependent variables are one)
margEffCondOne <- margEff( estResultBFGS, cond = TRUE, othDepVar = 1,
   algorithm = GenzBretz() )
print( margEffCondOne )

# conditional marginal effects with covariance matrix at sample mean
# (assuming that all other dependent variables are one)
margEffCondOneCov <- margEff( estResultBFGS, cond = TRUE, othDepVar = 1,
   data = as.data.frame( t( colMeans( xMat ) ) ), calcVCov = TRUE,
   algorithm = GenzBretz() )
print( margEffCondOneCov )
print( attr( margEffCondOneCov, "vcov" ) )
print( drop( attr( margEffCondOneCov, "vcov" ) ) )
summary( margEffCondOneCov )
# now with using argument 'atMean'
margEffCondOneCovA <- margEff( estResultBFGS, cond = TRUE, othDepVar = 1,
   atMean = TRUE, algorithm = GenzBretz(), calcVCov = TRUE )
all.equal( margEffCondOneCovA, margEffCondOneCov )
