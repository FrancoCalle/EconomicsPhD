#	Chebyshev-Approximator.R
#	-------------------------------------------------------------------------------------
#	Günter J. Hitsch, April 2014
#	Chebyshev approximation for arbitrary number of dimensions
#   Franco's Version in .jl


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


using Kronecker
using Distributions

mutable struct ChebyshevFeatures
    N_dim::Int
    N_degree::Int 
    N_nodes::Int
    lower_bound::Array 
    upper_bound::Array
    theta::Array
    X::Array
    TX::Array
    T2::Array
end


function initializeChebyshevApproximator(N_dim, N_degree, N_nodes, lower_bound, upper_bound)
   
    # Catch some basic input errors
	if N_nodes < N_degree + 1
		print("Error: Number of Chebyshev interpolation nodes (N_nodes) < N_degree + 1!\n")
    end
	
	if sum(lower_bound >= upper_bound) > 0
		print("Error: Rectangle (grid) upper_bound > lower_bound not true in all dimensions!\n")
    end

    # Calculate original [-1, 1] nodes
	nu_nodes = .- cos.(pi .* ( (2 .* Array(1:N_nodes) .- 1) ./ (2 * N_nodes)))

	# Create a table of all N_dim dimensional nodes (N_nodes^N_dim by N_dim matrix)
    list_nodes_nu = ()
    for z=1:N_dim
            list_nodes_nu = (list_nodes_nu..., nu_nodes)
    end
    X_temp = collect(Base.product(list_nodes_nu...))
    X_temp = reshape(X_temp, (length(nu_nodes) ^ N_dim))    
    
    # Placeholder
    X = zeros(length(nu_nodes) ^ N_dim, N_dim)
    for ii in 1:length(nu_nodes) ^ N_dim
        X[ii,:] = collect(X_temp[ii])'
    end
	X = X[:, end:-1:1] # Flip.

    # Adjust nodes to range of the interval in dimension k
	delta = upper_bound .- lower_bound
    X = 0.5 .* delta' .* (X .+ 1) .+ lower_bound'

    # TODO:::
	# Calculate Chebyshev polynomial values for all nodes (N_nodes by N_degree+1 matrix)
	Tnu = calculateChebyshevPolynomials(N_degree, nu_nodes)
	
	# Calculate TX, the regressor matrix of Chebyshev polynomials
	TX = copy(Tnu)
	if N_dim > 1
		for k = 2:N_dim 
            TX = TX ⊗ Tnu
        end
    end
	
	# Precompute the denominator "sum of squares" terms used in the Chebyshev coefficients formula
	K = (N_degree + 1)^N_dim
    T2 = sum(TX .* TX, dims=1)'

	# Initialize Chebyshev regression coefficients
	theta = zeros(K, 1)
	
	cheb = ChebyshevFeatures(N_dim, 
                                N_degree, 
                                N_nodes,
                                lower_bound, 
                                upper_bound,
                                theta, 
                                X, 
                                TX, 
                                T2)

    return cheb

end


function calculateChebyshevPolynomials(N_degree, x)
	
	Tx = zeros(length(x), N_degree + 1)

	# Recursively calculate polynomials
	Tx[:,1] .= 1
	if N_degree >= 1 
        Tx[:,2] .= x					# If polynomial degree >= 1
    end
	if N_degree >= 2	            # If polynomial degree >= 2
		for k in 2:N_degree
			Tx[:,k+1] .= 2 .* x .* Tx[:,k] .- Tx[:,k-1]
        end
    end

	return Tx
end


function calculateChebyshevCoefficients(f, cheb) 
	
    theta_hat = (cheb.TX' * f(cheb.X)) ./ cheb.T2

    setfield!(cheb, :theta, theta_hat)

	return cheb
end


function evaluateChebyshev(Z, cheb) 
	
	L = size(Z)[1]
    K = cheb.N_dim
		
	# Adjust points in Z to rectangle
    Z = 2 .* (Z .- cheb.lower_bound' ) ./ (cheb.upper_bound .- cheb.lower_bound)' .- 1
	
	# Polynomial values for each row in Z
	TZ = zeros(L, (cheb.N_degree + 1)^cheb.N_dim)

    for i = 1:L 	
		z = calculateChebyshevPolynomials(cheb.N_degree, Z[i,K])
		if K > 1
			for k =  1:K-1 
                z = calculateChebyshevPolynomials(cheb.N_degree, Z[i,K-1 - (k-1)]) ⊗ z
            end
        end
		TZ[i,:] = z
    end
	
	y_hat = TZ * cheb.theta
	
	return y_hat
end
