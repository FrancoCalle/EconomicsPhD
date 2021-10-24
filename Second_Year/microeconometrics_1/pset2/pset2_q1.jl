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

#Load packages ...

using Distributions
using LinearAlgebra
# using StatsBase
# using CategoricalArrays
using DataFrames
using Plots
using StatFiles
using Chain
using Tables
using CSV
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
        uv = rand(MvNormal(μ_uv, σ_uv), N)'
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


# Now multiple weak IV
#-------------------------------------
Pkg.add("FixedEffectModels")
using FixedEffectModels

include(joinpath(@__DIR__,"..", "..", "..","fc_toolkit.jl"))

function multiple_iv_mc_simulations(γ, ρ_uv=0.99, n=1000, m=10000)
    
    k = length(γ)
    μ_uv = [0, 0]
    σ_uv = Matrix{Float64}(I, 2, 2); σ_uv[1,2] = ρ_uv; σ_uv[2,1] = ρ_uv
    ivresults = Vector{Float64}()
    
    for mm in 1:m
        uv = rand(MvNormal(μ_uv, σ_uv), N)'
        Σ_z = Diagonal(ones(k)) .+ zeros(k,k)
        Z = rand(MvNormal(zeros(k), Σ_z), n)'
        D = Z * γ .+ uv[:,2]
        Y = uv[:,1]
        df = DataFrame(:Y => Y[:], :D => D[:], :Z1 => Z[:,1][:])
        β = coef(reg(df, @formula(Y ~ 0 + (D ~ 0 + Z1 ))))
        append!(ivresults, β)
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
        println("Used cholesky decomposition...")
        μ  = γ'*U / σ_v
    else
        println("Used plain square mat...")
        μ  = γ'*sqrt(Q_zz) / σ_v
    end

    approximation = (σ_u / σ_v) .* (μ*rrf')./(μ²./√k .+ μ*rfs')
    
    # Filter to values from -2 to 4, as in class notes
    return filter(x -> (x > -3 && x < 3), approximation)

end  


ρ = 0.9
k = 100;
γ = ones(k,1).*0.025;

density(multiple_iv_mc_simulations(γ,ρ), xlims = (-3,3), label="True Dist",linewidth=2.5)
density!(multivariate_weak_iv_approximation(γ,ρ), xlims = (-3,3), label="Weak IV Asymptotics",linewidth=2.5, linestyle=:dash)
# density!(multivariate_weak_iv_approximation(γ,ρ), xlims = (-3,3), label="Weak IV Asymptotics",linewidth=2.5, linestyle=:dash)
density!(multivariate_standard_approximation(γ), xlims = (-3,3), label="Standard Asymptotics",linewidth=2.5, linestyle=:dash)


# Numerical example:
approximation = multivariate_weak_iv_approximation(γ ,ρ , 1000, 100000);
approximation_cholesky = multivariate_weak_iv_approximation(γ ,ρ , 1000, 100000, true);

density(approximation, xlims = (-3,3), label="√Q_zz",linewidth=2.5, color=:red)
density!(approximation_cholesky, xlims = (-3,3), label="Cholesky Decomp", linewidth=2.5, color=:blue, linestyle=:dash)






