include("Gauss-Hermite-Quadrature.jl")


using Plots
# Settings ------------------------------------------------------------------------------
N_dim 		= 3						# Dimensionality = number of function arguments
N_nodes		= 9						# No. of quadrature nodes
L			= 10000000 				# Number of simulation draws

mu = [0.5, -1.0, 0.0]				# Multivariate normal distribution: mean

# Covariance matrix: creative positive definite matrix (a valid covariance matrix) from A
A = Array([1.0  0.2 -0.5; 
        0.2  1.6  0.9; 
        -0.5  0.9  0.2])

# Generate 
Cov = A' * A			# Covariance matrix of multivariate normal distribution
chol_Cov = cholesky(Cov).U        # Cholesky decomposition


# Functions to integrate ----------------------------------------------------------------
f(X) = sin.(2 .* X[:,1]) .* X[:,2].^2 .+ X[:,3]


# Initialize quadrature nodes and weights
gauss_hermite = GaussHermite(N_dim, 50)

X, weight = gauss_hermite

# Gauss-Hermite integral
M			= size(X,1)
X_adjusted	= sqrt(2) .* X * chol_Cov .+ repeat(mu',M, 1)  # Distrib is ∼ N(\mu, A), with nodes we are representing the support
y			= f(X_adjusted)
I_gauss_hermite	= (1/pi^(N_dim/2)) .* weight' * y

plot(y, weight)
xlims!((-50,50))

# Simulated (Monte Carlo) integral

X_draws		= rand(Normal(0,1), L, N_dim)  * chol_Cov .+ repeat(mu', L, 1)
y			= f(X_draws)
I_simulated	= mean(y)


print("Done")
print(I_gauss_hermite[1],'\n')
print(I_simulated,'\n')
print( abs(I_gauss_hermite[1] - I_simulated),'\n')
