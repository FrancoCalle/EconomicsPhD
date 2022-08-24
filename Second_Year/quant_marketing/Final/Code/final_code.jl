using Pkg
include("Gauss-Hermite-Quadrature.jl")
include("Chebyshev-Approximator.jl")

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
end


struct ModelParameters
    ρ::Float16
    α::Float16
    γ::Float16 
    β::Float16 
    σ_v::Float16 
    σ_0::Float16 
    μ_0::Float16
    ϑ::Array
end

struct ChebyshevParameters
    N_dim::Int 	        # Number of dimensions (variables of the model)
    N_degree::Int       # Degree / order of the Chebyshev polynomial
    N_nodes::Int        # Number of Chebyshev interpolation nodes
    lower_bound::Array	# Lowest value for each variable
    upper_bound::Array  # Highest value of each variable
end


dims = ModelDimensions(1000, 2, 5, 10)
mp   = ModelParameters(.9,      # ρ
                        .79,    # α
                        .5,     # γ
                        .995,   # β
                        .5,     # σ_v
                        .5,     # σ_0
                        1,      # μ_0
                        [.5 .5] # ϑ
                        )

function compute_random_utility(ξ_ij, mp = mp)

    u_ijt = mp.γ .- exp.(  -ρ .* ξ_ij  ) 

    return u_ijt

end

function update_beta(σ_prior, σ_v)

    β = σ_prior.^2 ./ (σ_prior.^2 .+ σ_v^2)

    return β

end

function update_quality_mean(μ_prior, D_jt, v_ij, β)

    μ_post = μ_prior + D_jt .* β .* (.-μ_prior .+ v_ij) 

    return μ_post
end

function update_quality_variance(σ_prior_0, σ_prior, D_jt)

    σ_post = sqrt.(1 ./ ( 1 ./ σ_prior_0.^2 .+  sum(D_jt, dims=2) ./ (σ_prior.^2) ))

    return σ_post
end


function generate_model_simulation(mp, dims)

    maxiter= 200

    # Unpack Parameters:
    α = mp.α
    β = mp.β
    σ_v = mp.σ_v
    σ_0 = mp.σ_0
    μ_0 = mp.μ_0
    ϑ = mp.ϑ

    # Unpack Dimensions:
    J = dims.J
    nStates = 50
    nPolynomials = dims.nPolynomials
    T = maxiter


    # Generate fake prices:
    P = rand(J,1) # rand(J,nStates) Assume just one price state for now ...

    # Generate logit shocks:
    ϵ_jt = rand(Gumbel(0,1), J, T)

    # Get nodes and weights from quadrature assuming only state is quality:
    q_nodes, q_weights = GaussHermite(1, nStates)

    # Generate arrays where we will store simulations: 
    
    # 1. Prior mean and variance:
    μ_prior = zeros(J, T) ; μ_prior[:,1] .= μ_0
    σ_prior = zeros(J, T) ; σ_prior[:,1] .= σ_0

    # 2. Utility:
    u_jt = zeros(J, T)
    
    u_prime_jt = zeros(J, T)

    # 3. Choice:
    D_jt = zeros(J, T)
    
    V_all = zeros(J, T)

    tt = 1
    V0 = 1
    V_prime = 0
    
    Vdiff = 10

    # COMPUTING VALUE FUNCTION: INCLUDE LEARNING STUFF WITHIN THE FIXED POINT
    while (Vdiff > 0.00000001) & (tt < maxiter)

        V0 = V_prime

        # Draw quality shock at time t:
        v_ij = rand(Normal(0, σ_v.^2), (J,1))

        # Compute today's quality experience: 
        ξ_ij = v_ij + ϑ'

        # Compute beliefs about tomorrow's quality experience:        
        ξ_prime_ij = q_nodes .* σ_prior[:,tt]' .+ μ_prior[:,tt]'

        #1. Generate today's Random Utility:
        u_jt[:,tt] = compute_random_utility(ξ_ij) .-  α .* P .+ ϵ_jt[:,tt]

        for jj = 1:J
            #2. Apply Chebyshev's approximation to tomorrow's random utility:
            # Initialize Chebyshev approximator
            cheb = initializeChebyshevApproximator(1, 5, 9, [minimum(ξ_prime_ij[:,jj])], [maximum(ξ_prime_ij[:,jj])])

            # Calculate Chebyshev regression coefficients to approximate f
            cheb = calculateChebyshevCoefficients(compute_random_utility, cheb)

            # Approximate random utility using cheby polynomials (First case, only)
            u_prime_jt[jj, tt] = dot(evaluateChebyshev(ξ_prime_ij[:,jj], cheb), q_weights) - α * P[jj]
        
        end

        V_all[:,tt] = u_jt[:,tt] + β .* u_prime_jt[:,tt]

        D_jt[:,tt] = V_all[:,tt] .== maximum(V_all[:,tt])

        V_prime = maximum(V_all[:,tt])

        
        # Update Beta:
        β_update = update_beta(σ_prior[:,tt], σ_v)

        # Update Prior Mean for next period: 
        μ_prior[:, tt+1] = update_quality_mean(μ_prior[:,tt], D_jt[:,tt], v_ij, β_update) 

        # Update Prior Variance for next period: 
        σ_prior[:, tt+1] = update_quality_variance(σ_prior[:,1], σ_prior[:,tt], D_jt[:,1:tt]) 

        Vdiff = abs(V0 - V_prime)

        print("Iteration number",tt, "; V_diff = ", Vdiff,"\n")
        
        tt += 1
    end

    T = tt
    plot(V_all[1,1:T])
    plot!(V_all[2,1:T])

    plot(μ_prior[1,1:T])

    return
end





cheb_hyperparams = ChebyshevParameters(3, 5, 9, [-5, -2, -3], [ 2,  4,  3])

N_params    = 1						# Dimensionality = number of function arguments
N_q_nodes	= 50				    # No. of quadrature nodes
L			= 10000000 				# Number of simulation draws
    

vijs = q_nodes .* σ_prior[:,1]' .+ μ_prior[:,1]'



compute_random_utility(vijs, μ_prior[:,tt], mp)