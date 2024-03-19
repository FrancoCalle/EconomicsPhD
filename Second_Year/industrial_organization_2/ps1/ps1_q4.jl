using Pkg

Pkg.update()
#Install ...
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add("Tables")
Pkg.add("CSV")
Pkg.add("Optim")
Pkg.instantiate()


#Load packages ...

using Distributions
using LinearAlgebra
using DataFrames
using StatFiles
# using Chain
using Tables
using CSV
using Optim
using Random
# using Plots
using Base.Threads
# using PyCall


# Excercise 2
function unpackVariables(dataset)
    
    s_jt = dataset.shares
    p = dataset.p
    x = dataset.x
    z = dataset[:,6:end]
    m = dataset.market
    j = dataset.choice
    M = maximum(m);
    J = maximum(j);

    return s_jt , p, x , z, m, j, M, J

end


function blpShares(δ_jt, m, X, Γ_νi)

    # nDrawsVi = 50  

    M = maximum(m)

    s_blp_jt = zeros(J*M)

    @threads for ii = 1:M  # obtain shares conditional on each market
        
        u_jn = δ_jt[m .== ii] .+ X[m .== ii,:]*Γ_νi[ii]
    
        pr_jn = exp.(u_jn)./(1 .+ sum(exp.(u_jn), dims=1))
        
        s_blp_jt[m .== ii]= mean(pr_jn, dims=2)
    
    end

   
    return s_blp_jt
end



function blpShares_with_demographics(δ_jt, # Mean Utilities
                                        m, # Markets
                                        p, # Prices
                                        x, # Quality 
                                        Π, # Parameters
                                        Γ_νi # Unobserved Heterogeneity
                                        )

    # nDrawsVi = 50  

    M = maximum(m)

    X = hcat(p, x)  # Add constant here ...

    X = hcat(ones(size(X,1)), X)

    s_blp_jt = zeros(J*M)

    @threads for ii = 1:M  # obtain shares conditional on each market
        
        u_jn = δ_jt[m .== ii] .+ X[m .== ii,:] * Π * D .+ X[m .== ii,:]*Γ_νi[ii]
    
        pr_jn = exp.(u_jn)./(1 .+ sum(exp.(u_jn), dims=1))
        
        s_blp_jt[m .== ii]= mean(pr_jn, dims=2)
    
    end

   
    return s_blp_jt
end


function fixed_point_function(δ_jt_lag, Γ_νi ; 
                                s_jt=s_jt, 
                                X)

    δ_jt = δ_jt_lag + log.(s_jt) - log.(blpShares(δ_jt_lag, m, X, Γ_νi))

    return δ_jt

end


function inner_loop(Γ, 
                    δ_jt_lag; # Changing
                    vi=vi, 
                    s_jt=s_jt, 
                    X,
                    M=M)

    tol = 1 # tolerance

    Γ_νi = [Γ*vi[ii] for ii = 1:M]  # Weight the nodes...
    
    iter = 0

    while tol > 10^-14

        δ_jt = fixed_point_function(δ_jt_lag, Γ_νi; s_jt, X)
        # δ_jt = δ_jt_lag + log.(s_jt) - log.(blpShares(δ_jt_lag, m, p, x, Γ_νi))

        tol = sum((δ_jt - δ_jt_lag).^2)

        iter += 1
        # print("Euclidean Distance: ", tol, "\n", "Iteration: ", iter, "\n")

        δ_jt_lag = δ_jt

    end

    δ_jt = δ_jt_lag

    return δ_jt

end


function compute_alpha(num, denom)

    alpha = -num/denom

    return alpha

end


function squarem_inner_loop(Γ, δ_jt_lag; # Changing
                    vi=vi, 
                    s_jt=s_jt, 
                    X = X,
                    M=M)

    tol = 1 # tolerance

    Γ_νi = [Γ*vi[ii] for ii = 1:M]  # Weight the nodes...

    iter = 0

    while tol > 10^-17

        # Get the residuals 
        δ_jt = fixed_point_function(δ_jt_lag, Γ_νi; s_jt, X)
        q1 = δ_jt .- δ_jt_lag
        δ_jt_prime = fixed_point_function(δ_jt, Γ_νi; s_jt, X)
        q2 = δ_jt_prime .- δ_jt

        # Form quadratic terms
        α_denom = (q2-q1)'*(q2-q1);
        α_num = q1'*(q2-q1);

        # Get the step-size
        alpha = compute_alpha(α_num,α_denom);
        δ_jt = δ_jt_lag + 2 * alpha * q1 + alpha.^2 * (q2-q1);

        # Fixed point iteration beyond the quadratic step
        δ_jt_prime = fixed_point_function(δ_jt, Γ_νi; s_jt, X);
    
        tol = sum((δ_jt_prime - δ_jt).^2)

        iter += 1
        # print("Euclidean Distance: ", tol, "\n", "Iteration: ", iter, "\n")

        δ_jt_lag = δ_jt_prime

    end

    δ_jt = δ_jt_lag

    return δ_jt

end



function get_marginal_cost(p, η_jkm)

    P = reshape(p, J, M)

    I_J= Matrix(I,J,J)
    mc_jm = zeros(J,M)
    for mm = 1:M 
        mc_jm[:,mm] = P[:,mm] .+ inv(I_J.*η_jkm[:,:,mm]) * P[:,mm]
    end

    return mc_jm

end


function gmm_objective(parameters; vi=vi, p = p, x = x, z = z)

    # Create X matrix Add constant here ...
    X = hcat(p, x)  
    X = hcat(ones(size(X,1)), X) # Add constant here 

    # Unpack parameters...
    α = parameters[1]
    β = parameters[2]

    Γ = zeros(3,3);
    Γ[1,1] = parameters[3];
    Γ[2,1] = parameters[4];
    Γ[3,1] = parameters[5];

    Γ[1,2] = Γ[2,1];
    Γ[2,2] = parameters[6];
    Γ[3,2] = parameters[7];
    
    Γ[1,3] = Γ[3,1];
    Γ[2,3] = Γ[3,2];
    Γ[3,3] = 0;

    α0 = parameters[8]

    # Compute contraction and get mean utility...
    δ_jt_lag = rand(M*J);
    # δ_jt = inner_loop(Γ, δ_jt_lag; vi=vi, s_jt=s_jt, X = X)
    δ_jt = squarem_inner_loop(Γ, δ_jt_lag; vi=vi, s_jt=s_jt, X = X)
    # Compute GMM objective Function
    Z = Array(z)
    Ω = inv(Z'*Z);
    
    ξ = δ_jt .- α0 - α * p - β * x

    obj = ξ'*Z*Ω*Z'ξ;

    return obj

end



#Compute Ealsticities:

function blpShares_nu(δ_jt, m, p, x, Γ_νi)

    nDrawsVi = size(Γ_νi[1],2)  

    M = maximum(m)

    X = hcat(p, x)

    s_blp_jtn = zeros(J*M, nDrawsVi)

    @threads for ii = 1:M  # obtain shares conditional on each market
        
        u_jn = δ_jt[m .== ii] .+ X[m .== ii,:]*Γ_νi[ii]
    
        pr_jn = exp.(u_jn)./(1 .+ sum(exp.(u_jn), dims=1))
                
        s_blp_jtn[m .== ii,:] = pr_jn
    end

   
    return s_blp_jtn
end


function get_elasticities(α, p, pr_jtn, s_jt; nDrawsVi=nDrawsVi)

    P = reshape(p, J, M)

    S = reshape(s_jt, J, M)

    η_jkm = zeros(J,J,M)


    for n = 1:nDrawsVi

        Pr_vi = reshape(pr_jtn[:,n], J, M)

        for mm = 1:M
            for jj = 1:J
                for kk = 1:J
                    if jj == kk

                        η_jkm[jj,kk,mm] +=  (-P[jj,mm]/S[jj,mm]) * (α .* Pr_vi[jj,mm] .* (1 - Pr_vi[jj,mm]))./nDrawsVi 

                    else

                        η_jkm[jj,kk,mm] += (P[kk,mm]/S[jj,mm]) * (α .* Pr_vi[jj,mm] .* Pr_vi[kk,mm])./nDrawsVi

                    end
                end
            end
        end
    end

    return η_jkm
end




#-----------------------------------------------

# Execute code:...
dataset = DataFrame(CSV.File("ps1_ex4.csv"));
# dataset = dataset[dataset.choice .==3,:]
# dataset.choice .= 1
nDrawsVi = 20
vi = [rand(MultivariateNormal([0,0,0], Matrix(I,3,3)), nDrawsVi) for ii = 1:M];

s_jt , p, x , z, m, j, M, J = unpackVariables(dataset);
K = 5
D = rand(size(z,1), K)

# Generate Random Mean Utilities ...


param_init = abs.(rand(3+5))

param_init = [
    -1
    1.5
   -2.06447019028676
    0.9916525538325098
    0.5416127398793276
    1.102057732887503
   -0.4288897218541497
   -2.441356865814114
]


# Run the optimization process:
result = optimize(gmm_objective, param_init, NelderMead(), 
                    Optim.Options(outer_iterations = 10000,
                                    iterations=2500,
                                    show_trace=true,
                                    show_every=100,
                                    g_tol = 1e-15,
                                    )
                                    )

params_hat = Optim.minimizer(result)

# Potential Candidate:
params_candidate1 = param_init





#-----------------------------------------------
# Obtain elasticities:

α = params_candidate1[1]

Γ = ones(3,3);
Γ[1,1] = params_candidate1[3];
Γ[2,1] = params_candidate1[4];
Γ[3,1] = params_candidate1[5];

Γ[1,2] = Γ[2,1];
Γ[2,2] = params_candidate1[6];
Γ[3,2] = params_candidate1[7];

Γ[1,3] = Γ[3,1];
Γ[2,3] = Γ[3,2];
Γ[3,3] = 0;

Γ_νi = [Γ*vi[ii] for ii = 1:M]  # Weight the nodes...

δ_jt = inner_loop(Γ, rand(M*J));

δ_jt_squarem = squarem_inner_loop(Γ, rand(M*J));

sum(abs.(δ_jt.- δ_jt_squarem))

pr_jtn = blpShares_nu(δ_jt, m, p, x, Γ_νi)

η_jkm = get_elasticities(α, p, pr_jtn, s_jt)

η_jk = zeros(J,J)
for mm = 1:M η_jk += η_jkm[:,:,mm]./M end # Average across markets:

round.(η_jk, digits= 6)



# Question 4, avg prices and quality...


s_jt , p, x , z, m, j, M, J

print("Prices:", mean(reshape(p,J,M), dims=2))
print("Quality:", mean(reshape(x,J,M), dims=2))
print("Shares", mean(reshape(s_jt,J,M), dims=2))



#-----------------------------------------------

using Pkg
Pkg.add("Optim")
using Optim

# Define a simple objective function, e.g., Rosenbrock's banana function
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2

# Run the optimizer
result = optimize(f, [0.0, 0.0], NelderMead(), Optim.Options(show_trace=true))

# Print the result
println(result)
