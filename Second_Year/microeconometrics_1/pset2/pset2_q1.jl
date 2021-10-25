# If you like what you see, don't forget to follow or(and) give me a star ⋆ in GitHub :)
using Pkg

#Install ...
Pkg.activate(".") 
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add(["DataFrames","DataFramesMeta","Chain"])
Pkg.add("Plots")
Pkg.add("CategoricalArrays")
Pkg.add("StatFiles")
Pkg.add("Tables")
Pkg.add("CSV")
Pkg.add("Optim")
Pkg.add("StatsPlots")
Pkg.instantiate()
Pkg.add("FixedEffectModels")

#Load packages ...

using FixedEffectModels
using Distributions
using LinearAlgebra
# using StatsBase
# using CategoricalArrays
using DataFrames
using Plots
using StatFiles
using Chain
using Tables
# using CSV
using Optim
using Random
using StatsPlots

include(joinpath(@__DIR__,"..", "..", "..","fc_toolkit.jl"))

# Replicate Monte Carlo simulation from class notes

# IV usual limit approximation:
#-------------------------------------

function iv_mc_simulations(γ ,ρ_uv=0.99, n=1000, m=10000)

    μ_uv = [0, 0]
    σ_uv = Matrix{Float64}(I, 2, 2); σ_uv[1,2] = ρ_uv; σ_uv[2,1] = ρ_uv
    ivresults = Vector{Float64}()

    for mm in 1:m
        uv = rand(MvNormal(μ_uv, σ_uv), n)'
        Z = rand(Normal(0,1),n)
        D = γ * Z + uv[:,2]
        Y = uv[:,1]
        β, _ = tsls_regression(Y, D, Z, nothing, false)
        append!(ivresults, β)
    end
    
    # Filter to values from -2 to 4, as in class notes
    return filter(x -> (x > -2 && x < 4), ivresults)
end


n = 1000
m = 10000
γ = 0.25
ρ_uv = .99

density(iv_mc_simulations(0.25), xlims = (-2.5,4.5), label="γ = 0.25",linewidth=2.5)
density!(iv_mc_simulations(0.1), xlims = (-2.5,4.5), label="γ = 0.1",linewidth=2.5)
density!(iv_mc_simulations(0.025), xlims = (-2.5,4.5), label="γ = 0.025",linewidth=2.5)
density!(iv_mc_simulations(0), xlims = (-2.5,4.5), label = "γ = 0",linewidth=2.5)
vline!([0],color=:gray, linestyle=:dash, label="",linewidth=2.5)
vline!([1],color=:gray, linestyle=:dash, label="",linewidth=2.5)



# Weak instrument asymptotics:
#-----------------------------

num_simdraws = 100000

function weak_iv_approx(γ_n, ρ_uv=0.99, num_simdraws = 100000, n = 1000)

    μ_uv = [0, 0]
    σ_uv = Matrix{Float64}(I, 2, 2); σ_uv[1,2] = ρ_uv; σ_uv[2,1] = ρ_uv

    rrf_rfs = rand(MvNormal(μ_uv, σ_uv), num_simdraws)'
    γ = γ_n*sqrt(n)
    approx = rrf_rfs[:,1]./(γ .+ rrf_rfs[:,2])

    return filter(x -> (x > -2 && x < 4), approx)

end


# Standard asymptotics:
#-----------------------------

function standard_approx(γ_n, num_simdraws = 100000, n = 1000)

    σ = 1/(sqrt(n)*γ_n)
    approx = rand(Normal(0, σ), num_simdraws)

    return filter(x -> (x > -2 && x < 4), approx)

end

density(iv_mc_simulations(0.025), xlims = (-2.5,4.5), label="True Dist",linewidth=2.5)
density!(weak_iv_approx(0.025), xlims = (-2.5,4.5), label="Weak IV Asymptotics",linewidth=2.5)
density!(standard_approx(0.025), xlims = (-2.5,4.5), label="Standard Asymptotics",linewidth=2.5)

γ_n = 0.2

#------------------------------------------------------------------------------------------------
# Part D: Now multiple weak IV
#------------------------------------------------------------------------------------------------
include(joinpath(@__DIR__,"..", "..", "..","fc_toolkit.jl"))


function multiple_iv_mc_simulations(γ, ρ_uv=0.99, n=1000, m=10000)
    
    k = length(γ)
    μ_uv = [0, 0]
    σ_uv = Matrix{Float64}(I, 2, 2); σ_uv[1,2] = ρ_uv; σ_uv[2,1] = ρ_uv
    ivresults = Vector{Float64}()
    
    for mm in 1:m
        uv = rand(MvNormal(μ_uv, σ_uv), n)'
        Σ_z = Diagonal(ones(k)) .+ zeros(k,k)
        Z = rand(MvNormal(zeros(k), Σ_z), n)'
        D = Z * γ .+ uv[:,2]
        Y = uv[:,1]
        df = DataFrame(:Y => Y[:], :D => D[:], :Z1 => Z[:,1][:])
        β = coef(reg(df, @formula(Y ~ 0 + (D ~ 0 + Z1 ))))
        append!(ivresults, β)
        """
        Not using this one because the previous tsls estimator is faster,
        but results are the same. If you want to see you can uncomment the
        following two lines and run the code again.
        β, _= tsls_regression(Y, D, Z[:,1], nothing, false) 
        append!(ivresults, β[1])
        """
    end
    
    # Filter to values from -2 to 4, as in class notes
    return filter(x -> (x > -3 && x < 3), ivresults)
end


k = 100
density(multiple_iv_mc_simulations(ones(k,1).*0.25, 0.5), xlims = (-2.5,4.5), label="γ = 0.25",linewidth=2.5)
density!(multiple_iv_mc_simulations(ones(k,1).*0.1, 0.5), xlims = (-2.5,4.5), label="γ = 0.1",linewidth=2.5)
density!(multiple_iv_mc_simulations(ones(k,1).*0.025, 0.5), xlims = (-2.5,4.5), label="γ = 0.025",linewidth=2.5)
density!(multiple_iv_mc_simulations(ones(k,1).*0), xlims = (-2.5,4.5), label = "γ = 0",linewidth=2.5)
vline!([0],color=:gray, linestyle=:dash, label="",linewidth=2.5)
vline!([1],color=:gray, linestyle=:dash, label="",linewidth=2.5)


# Asymptotic approximation:

# This will always be centered at 0 by definition
function multivariate_standard_approximation(γ_n, num_simdraws = 100000, n = 1000)

    k = length(γ_n)
    Σ = Diagonal((1 ./ γ_n')[:]).+ zeros((k,k))
    approx = rand(MvNormal(zeros(k), Σ), num_simdraws)'*γ_n

    # Filter to values from -2 to 4, as in class notes
    return filter(x -> (x > -3 && x < 3), approx)
end


function multivariate_weak_iv_approximation(γ_n ,ρ_uv=0.99, n=1000, num_simdraws = 100000, chol=false)

    σ_v = 1
    σ_u = 1
    γ = γ_n.*sqrt(n)
    k = length(γ_n)

    σ_uv = Matrix{Float64}(I, 2, 2); σ_uv[1,2] = ρ_uv; σ_uv[2,1] = ρ_uv
    Σ_z = Diagonal(ones(k)) .+ zeros(k,k)

    Z = rand(MvNormal(zeros(k), Σ_z), n)'
    Q_zz = (Z' * Z/n)
    

    Σ = kron(σ_uv,Q_zz)
    μ = zeros(size(Σ)[1])
    Φ = MvNormal(μ,Σ)
    rrf_rfs = rand(Φ, num_simdraws)'
   
    rrf = rrf_rfs[:,1:k]
    rfs =  rrf_rfs[:,k+1:end]

    μ² = γ'*Q_zz * γ / σ_v
    if chol == true
        U = Array(cholesky(Q_zz).U)
        μ  = γ'*U / σ_v
        # println("Used cholesky decomposition...")
    else
        μ  = γ'*sqrt(Q_zz) / σ_v
        # println("Used plain square mat...")
    end

    approximation = (σ_u / σ_v) .* (μ*rrf')./(μ²./√k .+ μ*rfs')
    
    # Filter to values from -2 to 4, as in class notes
    return filter(x -> (x > -3 && x < 3), approximation), μ²

end


# Numerical example:

K = 4
ρ = 0.9
γ = ones(K,1).*0.025;
true_distribution = multiple_iv_mc_simulations(γ,ρ)
approximation, _ = multivariate_weak_iv_approximation(γ ,ρ , 1000, 100000);
approximation_cholesky, _  = multivariate_weak_iv_approximation(γ ,ρ , 1000, 100000, true);

density(true_distribution, xlims = (-3,3), label="True Dist",linewidth=3, color=:blue)
density!(approximation_cholesky, xlims = (-3,3), label="Weak IV Approx. Cholesky",linewidth=2.5, linestyle=:dash)
density!(approximation, xlims = (-3,3), label="Weak IV Approx.",linewidth=3, linestyle=:dashdot)
density!(multivariate_standard_approximation(γ), xlims = (-3,3), label="Standard Asymptotics",linewidth=2.5, linestyle=:dash)
plot!(legend=:topleft)
savefig("Multivariate_IV_Approximation.pdf")


#------------------------------------------------------------
# Part E: Search μ^2 for different biases and K. 
#------------------------------------------------------------

function get_bias(δ, K=3, ρ=0.99)
    γ = ones(K,1).*δ;
    approximation, μ² = multivariate_weak_iv_approximation(γ , ρ , 1000, 100000);
    bias = mean(approximation)
    return [bias μ²[1]]
end

μ_grid = zeros(length(3:30), 4)
grid_δ = 0.001:0.0002:0.25

for K in 3:30
    println("Number of covariates: ", K)
    bias_list = get_bias.(grid_δ, K)
    bias_list = vcat(bias_list...)

    μ_grid[K-2, 1] = bias_list[argmin(abs.(bias_list[:,1] .- 0.05)),2]
    μ_grid[K-2, 2] = bias_list[argmin(abs.(bias_list[:,1] .- 0.1)),2]
    μ_grid[K-2, 3] = bias_list[argmin(abs.(bias_list[:,1] .- 0.2)),2]
    μ_grid[K-2, 4] = bias_list[argmin(abs.(bias_list[:,1] .- 0.3)),2]
end

CSV.write("q1E_table.csv",  Tables.table(μ_grid), writeheader=false)


#------------------------------------------------------------
# Part F: 
#------------------------------------------------------------

# Assume what stock and yogo said is true (TODO: derive the proof of it...)

c_value = zeros(length(3:30), 4)

for ii in 1:4
    for K in 3:30
        q = quantile(NoncentralChisq(K, μ_grid[K-2,ii].*K), .95)./K
        c_value[K-2, ii] = q
    end
end

CSV.write("c_value_table.csv",  Tables.table(c_value), writeheader=false)



plot(c_value[:,1], label=".05 Bias")
plot!(c_value[:,2], label=".1 Bias")
plot!(c_value[:,3], label=".2 Bias")
plot!(c_value[:,4], label=".3 Bias")
plot!(legend=:topleft)
savefig("CriticalValue_05.pdf")


plot(μ_grid[:,1], label=".05 Bias")
plot!(μ_grid[:,2], label=".1 Bias")
plot!(μ_grid[:,3], label=".2 Bias")
plot!(μ_grid[:,4], label=".3 Bias")
plot!(legend=:topleft)
savefig("Eigenvalue_bias.pdf")

#_--------------------------------------------------






















