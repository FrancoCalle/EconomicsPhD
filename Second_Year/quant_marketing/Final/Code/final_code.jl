using Pkg

# Pkg.instantiate()

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



struct ModelDimensions
    nObs::Int
    J::Int
    nStates::Int
    nPolynomials::Int
    # nChars::Int
end


struct ModelParameters
    ρ::Float16
    α::Float16
    γ::Float16 
    β::Float16 
    τ::Float16 
    σ_ν::Float16 
    σ_τ::Float16 
    σ_0::Float16 
    μ_0::Float16
    ϑ::Array
end

dims = ModelDimensions(1000, 2, 5, 10)
mp = ModelParameters(.9, .2, .5, .998, 1, 1, 1, 1, 5, [.5 .5])


function compute_random_utility(P, v_ij, μ_prior, mp)

    return u_ijt = mp.γ .- exp.(  -ρ .* (μ_prior + ν_ij )  ) .-  α .* P
end

function update_beta(σ_prior, σ_ν)

    β = σ_prior.^2 ./ (σ_prior.^2 .+ σ_ν^2)

    return β

end


function update_quality_mean(μ_prior, D_jt, ν_ij, β)

    μ_post = μ_prior + D_jt .* β .* (.-μ_prior .+ ν_ij) 

    return μ_post
end



function update_quality_variance(σ_prior, D_jt)

    σ_post = 1 ./ ( 1 ./ σ_prior .+  D_jt ./ σ_ν )

    return σ_post
end


function generate_model_simulation(mp, dims)

    # Unpack Parameters:
    ρ = mp.ρ
    α = mp.α
    γ = mp.γ
    β = mp.β
    τ = mp.τ
    σ_ν = mp.σ_ν
    # σ_τ = mp.σ_τ
    σ_0 = mp.σ_0
    μ_0 = mp.μ_0
    ϑ = mp.ϑ

    # Unpack Dimensions:
    J = dims.J
    nStates = dims.nStates
    nPolynomials = dims.nPolynomials
    T = 15

    # Generate fake prices:
    P = rand(J,nStates)

    # Generate logit shocks:
    ϵ_jt = rand(Gumbel(0,1), J, T)

    # Generate arrays where we will store simulations: 
    
    # 1. Prior mean and variance:
    μ_prior = zeros(J, T) ; μ_prior[:,1] .= μ_0
    σ_prior = zeros(J, T) ; σ_prior[:,1] .= σ_0

    # 2. Utility:
    u_jt = zeros(J, T)

    # 3. Probabilities:
    pr_jt = zeros(J+1, 1)
    
    # 4. Choice:
    D_jt = zeros(J, T)
    

    for tt = 1:T

        # Draw quality shock at time t:
        ν_ij = rand(Normal(0, σ_ν.^2), (J,1))

        # Generate Random Utility:
        u_jt[:,tt] = compute_random_utility(P, ν_ij, μ_prior[:,tt], mp) .+ ϵ_jt[:,tt]

        # Generate Probability and obtain D_jt:
        pr_jt[1:2] = exp.(u_jt[:,tt])./ (1 + sum(exp.(u_jt[:,tt])))
        pr_jt[3] = 1 - sum(pr_jt[1:2])

        # Get choice from probs:
        D_jt[:,tt] = (pr_jt .== maximum(pr_jt))[1:2]

        if tt != T

            # Update Beta:
            β = update_beta(σ_prior[:,tt], σ_ν)

            # Update Prior Mean for next period: 
            μ_prior[:, tt+1] = update_quality_mean(μ_prior[:,tt], D_jt[:,tt], ν_ij, β) 

            # Update Prior Variance for next period: 
            σ_prior[:, tt+1] = update_quality_variance(σ_prior[:, tt], D_jt[:,tt]) 

        end

    end


    return
end



function value_function_iteration()

    return
end


