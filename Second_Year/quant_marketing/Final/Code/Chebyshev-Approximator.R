#	Chebyshev-Approximator.R
#	-------------------------------------------------------------------------------------
#	Günter J. Hitsch, April 2014
#	Chebyshev approximation for arbitrary number of dimensions



#	initializeChebyshevApproximator
#	-------------------------------------------------------------------------------------
#	Inputs:
#	N_dim		... dimensionality of approximation problem
#	N_degree	... degree of Chebyshev polynomial
#	N_nodes		... number of Chebyshev interpolation nodes
#	lower_bound	... lower bounds of rectangle
#	upper_bound	... upper bounds of rectangle
#
#	Output: list "cheb"
#	theta		... Chebyshev coefficients for current approximation
#	X			... matrix of all nodes
#	TX			... regressor matrix of Chebyshev polynomials
#	T2			... pre-computed Chebyshev sum-of-squares

initializeChebyshevApproximator = function(N_dim, N_degree, N_nodes, lower_bound, upper_bound) {
		
	# Catch some basic input errors
	if (N_nodes < N_degree + 1)
		stop("Number of Chebyshev interpolation nodes (N_nodes) < N_degree + 1!\n")
	
	if (sum(lower_bound >= upper_bound) > 0)
		stop("Rectangle (grid) upper_bound > lower_bound not true in all dimensions!\n")
	
	
	# Calculate original [-1, 1] nodes
	nu_nodes = -cos(pi*((2*(1:N_nodes) - 1)/(2*N_nodes)))
	
	# Create a table of all N_dim dimensional nodes (N_nodes^N_dim by N_dim matrix)
	xi_nodes = list()
	for (k in 1:N_dim) xi_nodes[[k]] = nu_nodes
	X = as.matrix(expand.grid(xi_nodes))
	X = X[, N_dim:1]								# Flip: Represents order more intuitively
	
	# Adjust nodes to range of the interval in dimension k
	delta = upper_bound - lower_bound
	for (k in 1:N_dim)
		X[,k] = 0.5*delta[k]*(X[,k] + 1) + lower_bound[k]

	# Calculate Chebyshev polynomial values for all nodes (N_nodes by N_degree+1 matrix)
	Tnu = calculateChebyshevPolynomials(N_degree, nu_nodes)
	
	# Calculate TX, the regressor matrix of Chebyshev polynomials
	TX = Tnu
	if (N_dim > 1) {
		for (k in 2:N_dim) TX = TX %x% Tnu
	}
	
	# Precompute the denominator "sum of squares" terms used in the Chebyshev coefficients formula
	K = (N_degree + 1)^N_dim
	T2 = matrix(0, nrow = K, ncol = 1)
	for (k in 1:K)
		T2[k] = t(TX[,k]) %*% TX[,k]
	
	# Initialize Chebyshev regression coefficients
	theta = matrix(0, nrow = K, ncol = 1)
	
	cheb = list(N_dim = N_dim, N_degree = N_degree, N_nodes = N_nodes,
				lower_bound = lower_bound, upper_bound = upper_bound,
				theta = theta, X = X, TX = TX, T2 = T2)
		
	return(cheb)
}


#	calculateChebyshevPolynomials
#	-------------------------------------------------------------------------------------
#	Calculates the Chebyshev polynomials up to degree N_degree for each point in the
#	the vector x

calculateChebyshevPolynomials = function(N_degree, x) {
	
	Tx = matrix(0, nrow = length(x), ncol = N_degree + 1)

	# Recursively calculate polynomials
	Tx[,1] = 1
	if (N_degree >= 1) Tx[,2] = x					# If polynomial degree >= 1
	if (N_degree >= 2)	{							# If polynomial degree >= 2
		for (k in 2:N_degree)
			Tx[,k+1] = 2*x*Tx[,k] - Tx[,k-1]
	}

	return(Tx)
}


#	calculateChebyshevCoefficients
#	-------------------------------------------------------------------------------------
#	Approximation of function f.  Uses regression to update the Chebyshev coefficients

calculateChebyshevCoefficients = function(f, cheb) {
	
	cheb$theta = (t(cheb$TX) %*% f(cheb$X)) / cheb$T2
	
	return(cheb)
}


#	evaluateChebyshev
#	-------------------------------------------------------------------------------------
#	Predict approximated function values for each row in Z

evaluateChebyshev = function(Z, cheb) {
	
	L = dim(Z)[1]
	K = dim(Z)[2]
		
	# Adjust points in Z to rectangle
	for (k in 1:K)
		Z[,k] = 2*(Z[,k]-cheb$lower_bound[k])/(cheb$upper_bound[k] - cheb$lower_bound[k]) - 1
	
	# Polynomial values for each row in Z
	TZ = matrix(0, nrow = L, ncol = (cheb$N_degree + 1)^cheb$N_dim)
	for (i in 1:L) {	
		z = calculateChebyshevPolynomials(cheb$N_degree, Z[i,K])
		if (K > 1) {
			for (k in (K-1):1) z = calculateChebyshevPolynomials(cheb$N_degree, Z[i,k]) %x% z
		}
		TZ[i,] = z
	}
	
	y_hat = TZ %*% cheb$theta
	
	return(y_hat)
}

