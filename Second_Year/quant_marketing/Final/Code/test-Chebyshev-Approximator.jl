#test-Chebyshev-Approximator.jl (Works fine...) 
include("Chebyshev-Approximator.jl")


N_dim 		= 3
N_degree	= 5
N_nodes		= 9
lower_bound = [-5, -2, -3]		# Rectangle bounds
upper_bound = [ 2,  4,  3]
L			= 10000 			# Number of points to compare true function and approximation

f(X) = X[:,1] .* log.(5 .+ X[:,2]) .* X[:,3] 

# More difficult:
f(X) = X[:,1].^2 .* cos.( X[:,2] ) .* exp.(X[:,3])

Z = zeros(L,N_dim)

for k = 1:N_dim
	Z[:,k] = rand(Uniform(lower_bound[k], upper_bound[k]), L)
end

# Initialize Chebyshev approximator
cheb = initializeChebyshevApproximator(N_dim, N_degree, N_nodes, lower_bound, upper_bound)

# Calculate Chebyshev regression coefficients to approximate f
cheb = calculateChebyshevCoefficients(f, cheb)

# True and predicted values
y		= f(Z)
y_hat	= evaluateChebyshev(Z, cheb)

# Mean and max absolute difference between true and predicted values
diff		= abs.(y - y_hat)
mean_diff	= mean(diff)
max_diff	= maximum(diff)

# Print result
print(mean(diff))
print(maximum(diff))
