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
    X::Array #Market level characteristic
    Z::Array #Firm-Market characteristic
    Firm::Array # Firm index 
    Market::Array # Market index
end

struct EstimationData
    T::Int #Number of Markets
    K::Int #Number of Potential Entrants
    X::Array #Market level Characteristic
    Z::Array #Firm-Market characteristic
    I::Array #Indicator for entrants
    N::Array #Number of Entrants
    Firm::Array # Firm index
    Market::Array # Market index
end


function heterogeneous_profitability(z_it, η_t, ϵ_it; mp)

    γ = mp.ggamma; β = mp.β; α = mp.aalpha; δ = mp.δ; ρ = mp.rho

    ϕ_it = z_it * α + ρ * η_t + √(1-ρ^2) * ϵ_it
    
    return ϕ_it
end


function profit_function(x_t, z_it, ϕ_it, n_t; mp)
    
    γ = mp.ggamma; β = mp.β; α = mp.aalpha; δ = mp.δ; ρ = mp.rho

    v_t = γ + x_t * β - δ * log(n_t)

    π_it =  v_t + ϕ_it

    return π_it
end


function get_n_entrants(X, Z, Φ; mp , T, K)
    
    Π = zeros(K, T)
    
    for tt in 1:T

        n_t = 0
        
        x_t = X[tt] # Market level characteristic
                
        profitability_index = sortperm(Φ[:,tt], rev=true)

        for kk in 1:K
            
            nn = profitability_index[kk]

            z_it = Z[nn, tt]

            ϕ_it = Φ[nn, tt]

            Π[nn, tt] = profit_function(x_t, z_it, ϕ_it, n_t; mp)        

            n_t += 1 

        end
    end

    I = Π.>0
    
    N = sum(I, dims=1)

    return Π, I, N

end


function generate_nash_equilibrium(md, mp)

    # Unpack Data:
    X = md.X; Firms = md.Firm; Market = md.Market; T = md.T; K = md.K; Z = md.Z
    # Unpack Parameters:
    γ = mp.ggamma; β = mp.β; α = mp.aalpha; δ = mp.δ; ρ = mp.rho

    # Draw market specific shocks:
    Eta = rand(Normal(0,1), T)
    Epsilon = rand(Normal(0,1), T, K)

    Φ = heterogeneous_profitability.(Z[:], Eta[Market], Epsilon[:]; mp)
    Φ = reshape(Φ, K, T)

    Π, I, N = get_n_entrants(X, Z, Φ; mp, T, K)

    return Π, I, N

end


mp = ModelParameters(1,2,6,3,0.8)

T = 10; K = 30;
X = exp.(rand(Normal(0,1), T));
Z = rand(Normal(0,1), K , T);
Firm = vcat([collect(1:30) for ii = 1:T]...);
Market = vcat([Int.(ones(K)*ii) for ii = 1:T]...);


md = ModelData(T, K, X, Z, Firm, Market);
Π, I, N = generate_nash_equilibrium(md, mp)

estimationData = EstimationData(T, K, X, Z, I, N, Firm, Market)


function model_prediction(estimationData, mp; S=20)

    # Unpack Data:
    T = estimationData.T; 
    K = estimationData.K; 
    X = estimationData.X; 
    Z = estimationData.Z;
    I = estimationData.I;
    N = estimationData.N';
    Market= estimationData.Market

    # Unpack Parameters:
    γ = mp.ggamma; β = mp.β; α = mp.aalpha; δ = mp.δ; ρ = mp.rho

    # S draws:
    Eta_S = [rand(Normal(0,1), T, S) for ss = 1:S];
    Epsilon_S = [rand(Normal(0,1), T, K) for ss in 1:S];

    # S Simulations
    I_expected = zeros(K,T)
    N_expected = zeros(T)

    for ss = 1:S
        Φ = heterogeneous_profitability.(Z[:], Eta_S[ss][Market], Epsilon_S[ss][:]; mp)
        Φ = reshape(Φ, K, T)
        _, I_s, N_s = get_n_entrants(X, Z, Φ; mp, T, K)
        I_expected .+= I_s./S
        N_expected .+= N_s[:]./S
    end

    return I_expected, N_expected
    
end

