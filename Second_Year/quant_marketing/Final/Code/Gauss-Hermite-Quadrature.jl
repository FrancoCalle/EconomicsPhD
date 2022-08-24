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
using Kronecker
using Distributions

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

		z = 0
		z1 = -9
		its = 0
		p3 = 0.0
		pp = 0

        for i = 1:m 

            if i == 1
		        z = sqrt((2*n+1)) - 1.85575 * (2*n+1) ^ (-0.16667)
		    elseif i == 2
		        z = z - 1.14*(n^0.426)/z
		    elseif i == 3
		        z=1.86*z-0.86*x[1]
		    elseif i == 4
		        z=1.91*z-0.91*x[2]
		    else 
		        z=2.0*z-x[Int(i-2)]
            end

		    while (abs(z-z1) > EPS) | (its < MAXIT) 

				its += 1
		        p1 = PIM4
		        p2 = 0.0
		        
                for j = 1:n 
		            p3 = p2
		            p2 = p1
		            p1 = z * sqrt(2.0/j) * p2-sqrt(((j-1))/j) * p3
                end
		        
		        pp = sqrt(2*n)*p2
		        z1 = z
		        z = z1 - p1 / pp
				
            end
		    
		    x[Int(i)] = z
		    x[Int(n+1-i)] = -z
		    w[Int(i)] = 2.0 / (pp * pp)
		    w[Int(n+1-i)] = w[Int(i)]

        end
				
		return x, w
	end
	
	## Sub-function expand_mdim:
	function expand_mdim(x,y)
		
        N_x = size(x,1)
		N_y = size(y,1)
		z = 0
        
		for i = 1:N_x 
		    D = ones(N_y, size(x,2)) .* x[i,:]
		    E = hcat(D, y)
			if i==1
		    	z = copy(E)
			else
		    	z = vcat(z,E)
			end
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
	
	weight = cumprod(W, dims=2)
	weight = weight[:,dim]
	
	# gauss_hermite = list(X = X, weight = weight)

    return X, weight

end
