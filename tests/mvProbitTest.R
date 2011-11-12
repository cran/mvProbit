library( "mvProbit" )

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
sigma <- miscTools::symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

# all parameters in a vector
allCoef <- c( c( beta ), sigma[ lower.tri( sigma ) ] )

# generate dependent variables
yMatLin <- xMat %*% beta 
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
colnames( yMat ) <- paste( "y", 1:3, sep = "" )
# (yMatLin > 0 )== yMat

# unconditional expectations of dependent variables
yExp <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ) )
print( yExp )
yExpA <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = allCoef,
   data = as.data.frame( xMat ) )
all.equal( yExp, yExpA )
yExp2 <- pnorm( yMatLin )
all.equal( yExp, as.data.frame( yExp2 ) )


# conditional expectations of dependent variables
# (assuming that all other dependent variables are one)
yExpCond <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = GenzBretz() )
print( yExpCond )
yExpCondA <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = allCoef,
   data = as.data.frame( xMat ), cond = TRUE,
   algorithm = GenzBretz() )
all.equal( yExpCond, yExpCondA )
yExpCond2 <- matrix( NA, nrow = nObs, ncol = ncol( yMat ) )
for( i in 1:nObs ) {
   for( k in 1:ncol( yMat ) ) {
      set.seed( 123 )
      numerator <- pmvnorm( upper = yMatLin[ i, ], sigma = sigma )
      set.seed( 123 )
      denominator <- pmvnorm( upper = yMatLin[ i, -k ], sigma = sigma[ -k, -k ] )
      yExpCond2[ i, k ] <- numerator / denominator
   }
}
all.equal( yExpCond, as.data.frame( yExpCond2 ) )
# now with explicitly specifying the algorithm
yExpCond3 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = GenzBretz )
all.equal( yExpCond, yExpCond3 )
identical( yExpCond, yExpCond3 )
# now with integrals obtained by the Miwa algorithm
yExpCond4 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = Miwa )
all.equal( yExpCond, yExpCond4 )
# now with integrals obtained by the Miwa algorithm, less precise
yExpCond5 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = Miwa( steps = 32 ) )
all.equal( yExpCond4, yExpCond5 )
# now with integrals obtained by the TVPACK algorithm
yExpCond6 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = TVPACK )
all.equal( yExpCond, yExpCond6 )
# now with integrals obtained by the TVPACK algorithm, less precise
yExpCond7 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = TVPACK( abseps = 0.5 ) )
all.equal( yExpCond6, yExpCond7 )
# now with integrals obtained by the GHK algorithm
yExpCond8 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE )
all.equal( yExpCond, yExpCond8 )
# now with integrals obtained by the GHK algorithm, less precise
yExpCond9 <- mvProbitExp( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   nGHK = 100 )
all.equal( yExpCond8, yExpCond9 )


# conditional expectations of dependent variables
# (assuming that all other dependent variables are as observed)
yExpCondObs <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE, algorithm = GenzBretz() )
print( yExpCondObs )
yExpCondObsA <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE, algorithm = GenzBretz() )
all.equal( yExpCondObs, yExpCondObsA )
yExpCondObs2 <- matrix( NA, nrow = nObs, ncol = ncol( yMat ) )
for( i in 1:nObs ){
   for( k in 1:ncol( yMat ) ) {
      ySign <- 2 * yMat[ i, ] - 1
      ySign[ k ] <- 1
      yLinTmp <- yMatLin[ i, ] * ySign
      sigmaTmp <- diag( ySign ) %*% sigma %*% diag( ySign )
      set.seed( 123 )
      numerator <- pmvnorm( upper = yLinTmp, sigma = sigmaTmp )
      set.seed( 123 )
      denominator <- pmvnorm( upper = yLinTmp[ -k ], sigma = sigmaTmp[ -k, -k ] )
      yExpCondObs2[ i, k ] <- numerator / denominator
   }
}
all.equal( yExpCondObs, as.data.frame( yExpCondObs2 ) )
# now with explicitly specifying the algorithm
yExpCondObs3 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = GenzBretz )
all.equal( yExpCondObs, yExpCondObs3 )
identical( yExpCondObs, yExpCondObs3 )
# now with integrals obtained by the Miwa algorithm
yExpCondObs4 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = Miwa )
all.equal( yExpCondObs, yExpCondObs4 )
# now with integrals obtained by the Miwa algorithm, less precise
yExpCondObs5 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = Miwa( steps = 32 ) )
all.equal( yExpCondObs4, yExpCondObs5 )
# now with integrals obtained by the TVPACK algorithm
yExpCondObs6 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = TVPACK )
all.equal( yExpCondObs, yExpCondObs6 )
# now with integrals obtained by the TVPACK algorithm, less precise
yExpCondObs7 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = TVPACK( abseps = 0.5 ) )
all.equal( yExpCondObs6, yExpCondObs7 )
# now with integrals obtained by the GHK algorithm
yExpCondObs8 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE )
all.equal( yExpCondObs, yExpCondObs8 )
# now with integrals obtained by the GHK algorithm, less precise
yExpCondObs9 <- mvProbitExp( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, nGHK = 100 )
all.equal( yExpCondObs8, yExpCondObs9 )


# unconditional expectations of dependent variables by simulation
nSim <- 10000
ySim <- array( NA, c( nObs, ncol( yMat ), nSim ) )
for( s in 1:nSim ) {
   ySim[ , , s ] <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
}
yExpSim <- matrix( NA, nrow = nObs, ncol = ncol( yMat ) )
for( i in 1:nObs ) {
   yExpSim[ i, ] <- rowSums( ySim[ i, , ] ) / nSim
}
print( yExpSim )
print( yExpSim - as.matrix( yExp ) )

# for testing state of random number generator
rnorm( 4 )

# calculating log likelihood value(s)
logLikVal <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   algorithm = GenzBretz() )
print( logLikVal )
logLikValA <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ),
   algorithm = GenzBretz() )
all.equal( logLikVal, logLikValA )
# now with explicitly specifying the algorithm
logLikVal3 <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   algorithm = GenzBretz )
all.equal( logLikVal, logLikVal3 )
identical( logLikVal, logLikVal3 )
# now with integrals obtained by the Miwa algorithm
logLikVal4 <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   algorithm = Miwa )
all.equal( logLikVal, logLikVal4 )
# now with integrals obtained by the Miwa algorithm, less precise
logLikVal5 <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   algorithm = Miwa( steps = 32 ) )
all.equal( logLikVal4, logLikVal5 )
# now with integrals obtained by the TVPACK algorithm
logLikVal6 <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   algorithm = TVPACK )
all.equal( logLikVal, logLikVal6 )
# now with integrals obtained by the TVPACK algorithm, less precise
logLikVal7 <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   algorithm = TVPACK( abseps = 0.5 ) )
all.equal( logLikVal6, logLikVal7 )
# now with integrals obtained by the GHK algorithm
logLikVal8 <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ) )
all.equal( logLikVal, logLikVal8 )
# now with integrals obtained by the GHK algorithm, less precise
logLikVal9 <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   nGHK = 100 )
all.equal( logLikVal8, logLikVal9 )

# calculating log likelihood value(s) with one-sided gradients
logLikValGrad1 <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   oneSidedGrad = TRUE, algorithm = GenzBretz() )
print( logLikValGrad1 )
logLikValGrad1A <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ),
   oneSidedGrad = TRUE, algorithm = GenzBretz() )
all.equal( logLikValGrad1, logLikValGrad1A )

# calculating log likelihood value(s) with two-sided gradients
logLikValGrad <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   returnGrad = TRUE, algorithm = GenzBretz() )
print( logLikValGrad )
# now manually
llTmp <- function( coef ) {
   betaTmp <- coef[ 1:15 ]
   sigmaTmp <- diag( 3 )
   sigmaTmp[ lower.tri( sigmaTmp ) ] <- coef[ -(1:15) ]
   sigmaTmp[ upper.tri( sigmaTmp ) ] <- t( sigmaTmp )[ upper.tri( sigmaTmp ) ]
   result <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
      coef = betaTmp, sigma = sigmaTmp, 
      data = as.data.frame( cbind( xMat, yMat ) ), algorithm = GenzBretz() )
   return( result )
}
logLikValGrad2 <- numericGradient( llTmp, allCoef )
print( logLikValGrad2 )
attr( logLikValGrad1, "gradient" ) / logLikValGrad2 - 1
range( attr( logLikValGrad1, "gradient" ) / logLikValGrad2 - 1, na.rm = TRUE )
attr( logLikValGrad1, "gradient" ) - logLikValGrad2
range( attr( logLikValGrad1, "gradient" ) - logLikValGrad2 )
attr( logLikValGrad1, "gradient" ) / logLikValGrad2 - 1
range( attr( logLikValGrad, "gradient" ) / logLikValGrad2 - 1, na.rm = TRUE )
attr( logLikValGrad, "gradient" ) - logLikValGrad2
range( attr( logLikValGrad, "gradient" ) - logLikValGrad2 )

# for testing state of random number generator
rnorm( 4 )

# calculating marginal effects, unconditional
margEffUnc <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), vcov = diag( 18 ) )
print( margEffUnc )
print( attr( margEffUnc, "vcov" )[ 1:3, , ] )
print( drop( attr( margEffUnc, "vcov" )[ nObs, , ] ) )
summary( margEffUnc )
margEffUncA <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = allCoef,
   data = as.data.frame( xMat ), vcov = diag( 18 ) )
all.equal( margEffUnc, margEffUncA )
# now with explicitly specifying dummy variables
margEffUncD <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), vcov = diag( 18 ),
   dummyVar = c( "x1", "x2" ) )
all.equal( margEffUncD, margEffUnc )
# now with seemingly no dummy variables
margEffUncD0 <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), vcov = diag( 18 ),
   dummyVar = NULL )
summary( margEffUncD0 )
# now with seemingly only dummy variables
margEffUncDA <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), vcov = diag( 18 ),
   dummyVar = c( "x1", "x2", "x3", "x4" ) )
summary( margEffUncDA )
# now with mean values of the marginal effects
margEffUncM <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), vcov = diag( 18 ), 
   addMean = TRUE )
all.equal( margEffUnc, margEffUncM[ 1:nObs, ], check.attributes = FALSE )
print( margEffUncM[ nObs:(nObs+1), ] )
all.equal( attr( margEffUnc, "vcov" ), 
   attr( margEffUncM, "vcov" )[ 1:nObs, , ] )
print( attr( margEffUncM, "vcov" )[ nObs:(nObs+1), , ] )
print( drop( attr( margEffUncM, "vcov" )[ nObs+1, , ] ) )
all.equal( summary( margEffUnc )[ , ], 
   summary( margEffUncM )[ 1:( 12 * nObs ), ], check.attributes = FALSE )
printCoefmat( summary( margEffUncM )[ -( 1:( 12 * (nObs-1) ) ), ] )

# for testing state of random number generator
rnorm( 4 )

# calculating marginal effects, conditional
# (assuming that all other dependent variables are one)
margEffCond <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = GenzBretz() )
print( margEffCond )
margEffCondA <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = allCoef,
   data = as.data.frame( xMat ), cond = TRUE,
   algorithm = GenzBretz() )
all.equal( margEffCond, margEffCondA )
# now with dummy variables specified explicitly
margEffCondD <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   dummyVars = c( "x1", "x2" ), algorithm = GenzBretz() )
all.equal( margEffCondD, margEffCond )
# now with seemingly no dummy variables
margEffCondD0 <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   dummyVars = NULL, algorithm = GenzBretz() )
print( margEffCondD0 )
# now with seemingly only dummy variables
margEffCondDA <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   dummyVars = c( "x1", "x2", "x3", "x4" ), algorithm = GenzBretz() )
print( margEffCondDA )
# now with integrals obtained by the Miwa algorithm, reduced precision
margEffCond1 <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   algorithm = Miwa( steps = 32 ) )
print( margEffCond1 )
all.equal( margEffCond, margEffCond1 )
# now with integrals obtained by the GHK algorithm
margEffCond2 <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE )
print( margEffCond2 )
all.equal( margEffCond, margEffCond2 )
all.equal( margEffCond1, margEffCond2 )
# now with integrals obtained by the GHK algorithm, reduced precision
margEffCond3 <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat ), cond = TRUE,
   nGHK = 100 )
print( margEffCond3 )
all.equal( margEffCond, margEffCond3 )
all.equal( margEffCond2, margEffCond3 )
# now with variance covariance matrix and Jacobian
margEffCondV <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat )[ c(1,5,10), ], cond = TRUE, 
   vcov = diag( 18 ), returnJacobian = TRUE, algorithm = GenzBretz() )
print( attr( margEffCondV, "vcov" ) )
print( drop( attr( margEffCondV, "vcov" )[ 1, , ] ) )
print( attr( margEffCondV, "jacobian" ) )
print( drop( attr( margEffCondV, "jacobian" )[ 1, , ] ) )
summary( margEffCondV )
margEffCondVA <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = allCoef,
   data = as.data.frame( xMat )[ c(1,5,10), ], cond = TRUE, 
   vcov = diag( 18 ), returnJacobian = TRUE, algorithm = GenzBretz() )
all.equal( margEffCondV, margEffCondVA )
# now with Jacobian but without variance covariance matrix
margEffCondJac <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = allCoef,
   data = as.data.frame( xMat )[ c(1,5,10), ], cond = TRUE, 
   returnJacobian = TRUE, algorithm = GenzBretz() )
all.equal( attr( margEffCondJac, "jacobian" ), attr( margEffCondV, "jacobian" ) )
# now with mean values of marginal effects
margEffCondM <- mvProbitMargEff( ~ x1 + x2 + x3 + x4, coef = c( beta ), 
   sigma = sigma, data = as.data.frame( xMat )[ c(1,5,10), ], cond = TRUE, 
   vcov = diag( 18 ), addMean = TRUE, returnJacobian = TRUE,
   algorithm = GenzBretz() )
all.equal( margEffCondV, margEffCondM[ 1:3, ], check.attributes = FALSE )
print( margEffCondM )
all.equal( attr( margEffCondV, "vcov" ), 
   attr( margEffCondM, "vcov" )[ 1:3, , ] )
print( attr( margEffCondM, "vcov" ) )
print( drop( attr( margEffCondM, "vcov" )[ 4, , ] ) )
all.equal( attr( margEffCondV, "jacobian" ), 
   attr( margEffCondM, "jacobian" )[ 1:3, , ] )
print( attr( margEffCondM, "jacobian" ) )
print( drop( attr( margEffCondM, "jacobian" )[ 4, , ] ) )
all.equal( summary( margEffCondV )[ , ], summary( margEffCondM )[ 1:36, ],
   check.attributes = FALSE ) 
summary( margEffCondM )[ -( 1:24 ), ]

# for testing state of random number generator
rnorm( 4 )

# calculating marginal effects, conditional
# (assuming that all other dependent variables are as observed)
margEffCondObs <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE, algorithm = GenzBretz() )
print( margEffCondObs )
margEffCondObsA <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) ), 
   cond = TRUE, algorithm = GenzBretz() )
all.equal( margEffCondObs, margEffCondObsA )
# now with integrals obtained by the Miwa algorithm, reduced precision
margEffCondObs1 <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE, algorithm = Miwa( steps = 32 ) )
print( margEffCondObs1 )
all.equal( margEffCondObs, margEffCondObs1 )
# now with integrals obtained by the GHK algorithm
margEffCondObs2 <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, data = as.data.frame( cbind( xMat, yMat ) ),
   cond = TRUE )
print( margEffCondObs2 )
all.equal( margEffCondObs, margEffCondObs2 )
all.equal( margEffCondObs1, margEffCondObs2 )
# now with variance covariance matrix and Jacobian
margEffCondObsV <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) )[ c(1,5,10), ], 
   cond = TRUE, vcov = diag( 18 ), returnJacobian = TRUE,
   algorithm = GenzBretz() )
print( attr( margEffCondObsV, "vcov" ) )
print( drop( attr( margEffCondObsV, "vcov" )[ 1, , ] ) )
print( attr( margEffCondObsV, "jacobian" ) )
print( drop( attr( margEffCondObsV, "jacobian" )[ 1, , ] ) )
summary( margEffCondObsV )
margEffCondObsVA <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = allCoef, data = as.data.frame( cbind( xMat, yMat ) )[ c(1,5,10), ], 
   cond = TRUE, vcov = diag( 18 ), returnJacobian = TRUE,
   algorithm = GenzBretz() )
all.equal( margEffCondObs, margEffCondObsA )
# now with Jacobian but without variance covariance matrix
margEffCondObsJac <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) )[ c(1,5,10), ], 
   cond = TRUE, returnJacobian = TRUE, algorithm = GenzBretz() )
all.equal( attr( margEffCondObsJac, "jacobian" ), 
   attr( margEffCondObsV, "jacobian" ) )
# now with mean values of marginal effects
margEffCondObsM <- mvProbitMargEff( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4, 
   coef = c( beta ), sigma = sigma, 
   data = as.data.frame( cbind( xMat, yMat ) )[ c(1,5,10), ], 
   cond = TRUE, vcov = diag( 18 ), addMean = TRUE, returnJacobian = TRUE,
   algorithm = GenzBretz() )
all.equal( margEffCondObsV, margEffCondObsM[ 1:3, ], check.attributes = FALSE )
print( margEffCondObsM )
all.equal( attr( margEffCondObsV, "vcov" ), 
   attr( margEffCondObsM, "vcov" )[ 1:3, , ] )
print( attr( margEffCondObsM, "vcov" ) )
print( drop( attr( margEffCondObsM, "vcov" )[ 4, , ] ) )
all.equal( attr( margEffCondObsV, "jacobian" ), 
   attr( margEffCondObsM, "jacobian" )[ 1:3, , ] )
print( attr( margEffCondObsM, "jacobian" ) )
print( drop( attr( margEffCondObsM, "jacobian" )[ 4, , ] ) )
all.equal( summary( margEffCondObsV )[ , ], summary( margEffCondObsM )[ 1:36, ],
   check.attributes = FALSE )
summary( margEffCondObsM )[ -( 1:24 ), ]

# for testing state of random number generator
rnorm( 4 )

