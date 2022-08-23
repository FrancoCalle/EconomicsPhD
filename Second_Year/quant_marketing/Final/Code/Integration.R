#	Integration.R
#	-------------------------------------------------------------------------------------
#	Günter J. Hitsch, May 2014
#	Uses Gauss-Hermite quadrature to calculate an integral and compares with Monte Carlo
#	integration method

rm(list = ls())
source("Gauss-Hermite-Quadrature.R")



# Settings ------------------------------------------------------------------------------
N_dim 		= 3						# Dimensionality = number of function arguments
N_nodes		= 9						# No. of quadrature nodes
L			= 10000000 				# Number of simulation draws

mu = c(0.5, -1.0, 0.0)				# Multivariate normal distribution: mean

# Covariance matrix: creative positive definite matrix (a valid covariance matrix) from A
A	= matrix( c(  1.0,  0.2, -0.5,
				  0.2,  1.6,  0.9,
				 -0.5,  0.9,  0.2),  nrow=3, ncol=3, byrow=TRUE)

Cov = t(A) %*% A					# Covariance matrix of multivariate normal distribution
chol_Cov = chol(Cov)				# Cholesky decomposition



# Functions to integrate ----------------------------------------------------------------
f = function(X) {
	x1 = X[,1]
	x2 = X[,2]
	x3 = X[,3]
	
	return(sin(2*x1)*x2^2 + x3)
}
	


# Initialize quadrature nodes and weights
gauss_hermite = GaussHermite(N_dim, N_nodes)


# Gauss-Hermite integral
M			= dim(gauss_hermite$X)[1]
X_adjusted	= sqrt(2) * gauss_hermite$X %*% chol_Cov + matrix(mu, nrow=M, ncol=N_dim, byrow=TRUE)
y			= f(X_adjusted)
I_gauss_hermite	= (1/pi^(N_dim/2)) * t(gauss_hermite$weight) %*% y


# Simulated (Monte Carlo) integral
set.seed(124)

X_draws		= matrix(rnorm(L*N_dim), nrow=L, ncol=N_dim)  %*% chol_Cov + matrix(mu, nrow=L, ncol=N_dim, byrow=TRUE)
y			= f(X_draws);
I_simulated	= mean(y);


comparison = c(I_gauss_hermite, I_simulated, abs(I_gauss_hermite-I_simulated))
print(comparison)
