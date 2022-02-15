using Pkg

Pkg.instantiate()

#Load packages ...

using Distributions
using LinearAlgebra
using DataFrames
using StatFiles
using Tables
using CSV
using Optim
using Random
using Plots

struct ModelParameters
    aalpha::Float16 # Firm-market specific characteristic
    β::Float16 # Elasticity
    δ::Float16 # Propensity over number of entrants
    ggamma::Float16 #Profit Intercept
    rho::Float16 # Error Autocorrelation
end

struct ModelData
    T::Int #Number of Markets
    K::Int #Number of Potential Entrants
    X::Array #Number of Markets
    Z::Array #Number of Potential Entrants
    Firm::Array # Firm-market specific characteristic
    Market::Array # Elasticity
end


function heterogeneous_profitability(z_it, η_it, ϵ_it; mp)

    γ = mp.ggamma; β = mp.β; α = mp.aalpha; δ = mp.δ; ρ = mp.rho

    ϕ_it = z_it * α + ρ * η_t + √(1-ρ^2) * ϵ_it
    
    return ϕ_it
end


function profit_function(x_t, z_it, ϵ_it, η_t; mp)
    
    γ = mp.ggamma; β = mp.β; α = mp.aalpha; δ = mp.δ; ρ = mp.rho

    ϕ_it = heterogeneous_profitability(z_it, ϵ_it, mp)

    v_t = γ + x_t * β - δ * log(n_t)

    π_it =  v_t + ϕ_it

    return π_it
end


function generate_nash_equilibrium(md, mp)

    # Unpack Data:
    X = md.X; Z = md.Z; Firms = md.Firm; Market = md.Market; T = md.T; K = md.K
    # Unpack Parameters:
    γ = mp.ggamma; β = mp.β; α = mp.aalpha; δ = mp.δ; ρ = mp.rho

    # Draw market specific shocks:
    Eta = rand(Normal(0,1), T)
    Epsilon = rand(Normal(0,1), T*K)

    Φ = heterogeneous_profitability.(Z, Eta[Market], Epsilon; mp)
    Π = zeros(T*K)

    for jj in 1:T*K

        tt = Market[jj]

        x_t = X[tt]
        η_t = Eta[tt]

        z_it = Z[jj]
        ϵ_it = Epsilon[jj]

        Π[jj] = profit_function(x_t, z_it, ϵ_it, η_t; mp)        

    end


    return n_star

end



mp = ModelParameters(1,2,6,3,0.8)

T = 10; K = 30;
X = exp.(rand(Normal(0,1), K));
Z = rand(Normal(0,1), K * T);
Firm = vcat([collect(1:30) for ii = 1:T]...);
Market = vcat([Int.(ones(K)*ii) for ii = 1:T]...);


md = ModelData(T, K, X, Z, Firm, Market);





