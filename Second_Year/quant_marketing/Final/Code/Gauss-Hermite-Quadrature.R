#	Gauss-Hermite-Quadrature.R
#	-------------------------------------------------------------------------------------
#	Günter J. Hitsch, May 2014
#	Gauss-Hermite quadrature 



#	GaussHermite
#	-------------------------------------------------------------------------------------
#	Inputs:
#	N_dim		... dimensionality of integration problem
#	N_nodes		... number of nodes
#
#	Output: list "gauss_hermite"
#	X			... matrix of all nodes
#	weights		... 

GaussHermite = function(N_dim, N_nodes) {

	dim = N_dim
	K	= N_nodes
	
	## Sub-function gauher_mdim
	## (Exact copy of gauher)
	## Returns result as a list of x and w
	gauher_mdim = function(n) {
		EPS = 3.0e-14
		PIM4 = 0.7511255444649425
		MAXIT = 10
		x = matrix(0,nrow=n,ncol=1)
		w = matrix(0,nrow=n,ncol=1)
		m=(n+1)/2
		for (i in 1:m) {
		    if (i == 1)
		        z=sqrt((2*n+1))-1.85575*(2*n+1)^(-0.16667)
		    else if (i == 2)
		        z = z - 1.14*(n^0.426)/z
		    else if (i == 3)
		        z=1.86*z-0.86*x[1]
		    else if (i == 4)
		        z=1.91*z-0.91*x[2]
		    else 
		        z=2.0*z-x[i-2]
	
		    for (its in 1:MAXIT) {
		        p1=PIM4
		        p2=0.0
		        for (j in 1:n) {
		            p3=p2
		            p2=p1
		            p1=z*sqrt(2.0/j)*p2-sqrt(((j-1))/j)*p3
		        }
		        
		        pp=sqrt(2*n)*p2
		        z1=z
		        z=z1-p1/pp
		        if (abs(z-z1) <= EPS) break 
		    }
		    
		    x[i] = z
		    x[n+1-i] = -z
		    w[i]=2.0/(pp*pp)
		    w[n+1-i]=w[i]
		}
		
		result = list()
		result$x = x
		result$w = w
		
		return(result)
	}
	
	## Sub-function expand_mdim
	expand_mdim = function(x,y) {
		N_x = nrow(x)
		N_y = nrow(y)
		z = c()
		for (i in 1:N_x) {
		    D = matrix(x[i,], nrow=N_y, ncol=ncol(x), byrow=T)
		    E = cbind(D, y)
		    z = rbind(z, E)
		}
		
		return(z)
	}
	
	## gauher_multidim(dim ... dimensions, K ... quadrature nodes)
	## Creates an array X of quadrature node/vectors, and a vector of
	## corresponding weights
	
	result = gauher_mdim(K)
	
	## Create a list of nodes and weights
	X = result$x
	W = result$w
	for (i in 2:dim) {
	    X = expand_mdim(result$x,X)
	    W = expand_mdim(result$w,W)
	}
	
	weight = t(apply(W, 1, cumprod))
	weight = matrix(weight[,dim], ncol=1)
	
	gauss_hermite = list(X = X, weight = weight)
	return(gauss_hermite)
}

