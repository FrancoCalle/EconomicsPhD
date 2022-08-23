#	Chebyshev-Approximator.R
#	-------------------------------------------------------------------------------------
#	Günter J. Hitsch, April 2014
#	Chebyshev approximation for arbitrary number of dimensions
#   Franco's Version:



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



reshape([ x * y * z  for x=nu_nodes, y=nu_nodes, z=nu_nodes ], (length(nu_nodes) ^ 3, 9) )

function initializeChebyshevApproximator(N_dim, N_degree, N_nodes, lower_bound, upper_bound)
   
    # Catch some basic input errors
	if N_nodes < N_degree + 1
		print("Error: Number of Chebyshev interpolation nodes (N_nodes) < N_degree + 1!\n")
	
	if sum(lower_bound >= upper_bound) > 0
		print("Error: Rectangle (grid) upper_bound > lower_bound not true in all dimensions!\n")

    # Calculate original [-1, 1] nodes
	nu_nodes = .- cos.(pi .* ( (2 .* Array(1:N_nodes) .- 1) ./ (2 * N_nodes)))

	# Create a table of all N_dim dimensional nodes (N_nodes^N_dim by N_dim matrix)
    X_temp = collect(Base.product(xi_nodes[:,1], xi_nodes[:,2], xi_nodes[:,3]))
    X_temp = reshape(X_temp, (length(nu_nodes) ^ N_dim))    
    
    # Placeholder
    X = zeros(length(nu_nodes) ^ N_dim, N_dim)
    for ii in 1:length(nu_nodes) ^ N_dim
        X[ii,:] = collect(X_temp[ii])'
    end
	X = X[:, end:-1:1] # Flip.

    # TODO:::
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

end


N_dim 		= 3
N_degree	= 5
N_nodes		= 9
lower_bound = (-5, -2, -3)		# Rectangle bounds
upper_bound = ( 2,  4,  3)
L			= 10000 			# Number of points to compare true function and approximation
