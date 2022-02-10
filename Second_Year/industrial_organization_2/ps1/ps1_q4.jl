using Pkg

#Install ...
# Pkg.add("Distributions")
# Pkg.add("StatsBase")
# Pkg.add("Plots")
# Pkg.add("StatFiles")
# Pkg.add("Tables")
# Pkg.add("CSV")
# Pkg.add("Optim")
Pkg.instantiate()

#Load packages ...

using Distributions
using LinearAlgebra
using DataFrames
using StatFiles
using Chain
using Tables
using CSV
using Optim
using Random
using Plots
using Base.Threads
using PyCall
# pyblp = pyimport("pyblp")

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


function blpShares(δ_jt, m, p, x, Γ_νi)

    # nDrawsVi = 50  

    M = maximum(m)

    X = hcat(p, x)

    s_blp_jt = zeros(J*M)

    @threads for ii = 1:M  # obtain shares conditional on each market
        
        u_jn = δ_jt[m .== ii] .+ X[m .== ii,:]*Γ_νi[ii]
    
        pr_jn = exp.(u_jn)./(1 .+ sum(exp.(u_jn), dims=1))
        
        s_blp_jt[m .== ii]= mean(pr_jn, dims=2)
    
    end

   
    return s_blp_jt
end

function inner_loop(Γ, δ_jt_lag; # Changing
                    vi=vi, 
                    s_jt=s_jt, 
                    p=p, 
                    x=x,
                    M=M)

    tol = 1 # tolerance

    Γ_νi = [Γ*vi[ii] for ii = 1:M]  # Weight the nodes...

    while tol > 10^-14

        δ_jt = δ_jt_lag + log.(s_jt) - log.(blpShares(δ_jt_lag, m, p, x, Γ_νi))

        tol = sum((δ_jt - δ_jt_lag).^2)

        # print("Euclidean Distance: ", tol, "\n")

        δ_jt_lag = δ_jt

    end

    δ_jt = δ_jt_lag

    return δ_jt

end



function get_elasticities(α, p, s_jt, vi)

    S = reshape(s_jt, J, M)
    P = reshape(p, J, M)

    η_jkm = zeros(J,J,M)

    for mm = 1:M
        for jj = 1:J
            for kk = 1:J
                if jj == kk
                    η_jkm[jj,kk,mm] = -α .* P[kk,mm] .* (1 - S[kk,mm]) #Check whether this is the right approach to ge elasticiities
                else
                    η_jkm[jj,kk,mm] = α .* P[kk,mm] .* S[kk,mm]
                end
            end
        end
    end

    return η_jkm
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

    α = parameters[1]
    β = parameters[2]

    Γ = zeros(2,2);
    Γ[1,1] = parameters[3];
    Γ[2,1] = parameters[4];
    Γ[2,2] = parameters[5];

    # Compute contraction and get mean utility...
    δ_jt_lag = rand(M*J);
    δ_jt = inner_loop(Γ, δ_jt_lag)

    # Compute GMM objective Function
    Z = Array(z)
    Ω = inv(Z'*Z);
    
    ξ = δ_jt - α * p - β * x

    obj = ξ'*Z*Ω*Z'ξ;

    return obj

end



# Execute code:...
dataset = DataFrame(CSV.File("ps1_ex4.csv"));
s_jt , p, x , z, m, j, M, J = unpackVariables(dataset);

# Generate Random Mean Utilities ...
nDrawsVi = 20
param_init = abs.(rand(5))
param_init[1] = -param_init[1]
vi = [rand(MultivariateNormal([0,0], Matrix(I,2,2)), nDrawsVi) for ii = 1:M];

# gmm_objective(param_init)
result = optimize(gmm_objective, param_init, NelderMead(), Optim.Options(outer_iterations = 1500,
                    iterations=10000,
                    show_trace=true,
                    show_every=100))

params_hat = Optim.minimizer(result)

# Potential Candidate:
params_candidate1 = [
    -0.18472568371716075
    2.1385463127358983
    -2.0019003342751214
    -0.04381409409480813
    0.14376694750798633
]

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


# Obtain elasticities:

α = -params_candidate1[1]

Γ = zeros(2,2);
Γ[1,1] = params_candidate1[3];
Γ[2,1] = params_candidate1[4];
Γ[2,2] = params_candidate1[5];
Γ_νi = [Γ*vi[ii] for ii = 1:M]  # Weight the nodes...

δ_jt = inner_loop(Γ, rand(M*J))

pr_jtn= blpShares_nu(δ_jt, m, p, x, Γ_νi)

η_jkm = get_elasticities(α, p, pr_jtn, s_jt)

η_jk = zeros(J,J)
for mm = 1:M η_jk += η_jkm[:,:,mm]./M end # Average across markets:

round.(η_jk, digits= 6)



# Question 4, avg prices and quality...


s_jt , p, x , z, m, j, M, J

print("Prices:", mean(reshape(p,J,M), dims=2))
print("Quality:", mean(reshape(x,J,M), dims=2))
print("Shares", mean(reshape(s_jt,J,M), dims=2))


