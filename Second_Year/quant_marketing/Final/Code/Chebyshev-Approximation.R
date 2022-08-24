#	Chebyshev-Approximation.R
#	-------------------------------------------------------------------------------------
#	Günter J. Hitsch, April 2014
#	Uses Chebyshev approximation to approximate some function f

rm(list = ls())
source("Chebyshev-Approximator.R")


# Functions to approximate

# Polynomial
#f = function(X) {
#	X1 = X[,1]
#	X2 = X[,2]
#	X3 = X[,3]
#	return( X1*X3^3 + X2*X3 + X1^2*X2*X3^2  )
#}

# Medium difficult
f = function(X) {
	return( X[,1]*log(5 + X[,2])*X[,3] )
}

# More difficult
#f = function(X) {
#	return( X[,1]^2*cos(X[,2])*exp(X[,3]) )
#}


# Settings
N_dim 		= 3
N_degree	= 5
N_nodes		= 9
lower_bound = c(-5, -2, -3)		# Rectangle bounds
upper_bound = c( 2,  4,  3)
L			= 10000 			# Number of points to compare true function and approximation


# Initialize Chebyshev approximator
cheb = initializeChebyshevApproximator(N_dim, N_degree, N_nodes, lower_bound, upper_bound)

# Calculate Chebyshev regression coefficients to approximate f
cheb = calculateChebyshevCoefficients(f, cheb)

# Matrix for comparison
set.seed(124)
Z = matrix(0, nrow = L, ncol = N_dim)
for (k in 1:N_dim)
	Z[,k] = runif(L, min = lower_bound[k], max = upper_bound[k])

# True and predicted values
y		= f(Z)
y_hat	= evaluateChebyshev(Z, cheb)

# Mean and max absolute difference between true and predicted values
diff		= abs(y - y_hat)
mean_diff	= mean(diff)
max_diff	= max(diff)

#print(cbind(y, y_hat))

print(mean(diff))
print(max(diff))
