#	Gauss-Hermite-Quadrature.R
#	-------------------------------------------------------------------------------------
#	Günter J. Hitsch, May 2014
#	Gauss-Hermite quadrature
#   Franco Calle's version in .jl



#	GaussHermite
#	-------------------------------------------------------------------------------------
#	Inputs:
#	N_dim		... dimensionality of integration problem
#	N_nodes		... number of nodes
#
#	Output: list "gauss_hermite"
#	X			... matrix of all nodes
#	weights		... 

using LinearAlgebra

function GaussHermite(N_dim, N_nodes) 

	dim = N_dim
	K	= N_nodes
	
	## Sub-function gauher_mdim
	## (Exact copy of gauher)
	## Returns result as a list of x and w
    function gauher_mdim(n) 

        EPS = 3.0e-14
		PIM4 = 0.7511255444649425
		MAXIT = 10

        x = zeros(n,1)
		w = zeros(n,1)
		m = (n+1)/2

        for i = 1:m 

            if (i == 1)
		        z = sqrt((2*n+1)) - 1.85575 * (2*n+1) ^ (-0.16667)
		    elseif (i == 2)
		        z = z - 1.14*(n^0.426)/z
		    elseif (i == 3)
		        z=1.86*z-0.86*x[1]
		    elseif (i == 4)
		        z=1.91*z-0.91*x[2]
		    else 
		        z=2.0*z-x[i-2]
            end
	
		    for its = 1:MAXIT

		        p1=PIM4
		        p2=0.0
		        
                for j = 1:n 
		            p3=p2
		            p2=p1
		            p1=z * sqrt(2.0/j) * p2-sqrt(((j-1))/j) * p3
                end
		        
		        pp = sqrt(2*n)*p2
		        z1 = z
		        z = z1-p1/pp
		        if abs(z-z1) <= EPS break end
            end
		    
		    x[i] = z
		    x[n+1-i] = -z
		    w[i]=2.0 / (pp * pp)
		    w[n+1-i]=w[i]

        end
				
		return x, w
	end
	
	## Sub-function expand_mdim:
	function expand_mdim(x,y)
		
        N_x = nrow(x)
		N_y = nrow(y)
		z = c()
		
        for i = 1:N_x 
		    D = matrix(x[i,:], nrow=N_y, ncol=ncol(x), byrow=T)
		    E = cbind(D, y)
		    z = rbind(z, E)
        end
		
		return z
	end
	
	## gauher_multidim(dim ... dimensions, K ... quadrature nodes)
	## Creates an array X of quadrature node/vectors, and a vector of corresponding weights
	x, w = gauher_mdim(K)
	
	## Create a list of nodes and weights
    X = copy(x)
    W = copy(w)

    for i = 2:dim 
	    X = expand_mdim(x, X)
	    W = expand_mdim(w, W)
    end
	
	weight = apply(W, 1, cumprod)'
	weight = matrix(weight[:,dim], ncol=1)
	
	gauss_hermite = list(X = X, weight = weight)

    return gauss_hermite

end






# source("Gauss-Hermite-Quadrature.R")



# Settings ------------------------------------------------------------------------------
N_dim 		= 3						# Dimensionality = number of function arguments
N_nodes		= 9						# No. of quadrature nodes
L			= 10000000 				# Number of simulation draws

mu = [0.5, -1.0, 0.0]				# Multivariate normal distribution: mean

# Covariance matrix: creative positive definite matrix (a valid covariance matrix) from A
A = Array([1.0  0.2 -0.5; 0.2  1.6  0.9; -0.5  0.9  0.2])

# Generate 
Cov = A' * A					# Covariance matrix of multivariate normal distribution
chol_Cov = cholesky(Cov)        # Cholesky decomposition



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
