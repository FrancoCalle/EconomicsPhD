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
    T::Int
end


mutable struct ModelParameters
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
end

struct QuadratureParameters
    N_dim::Int 	        # Number of dimensions (variables of the model)
    N_States::Int       # Degree / order of the Chebyshev polynomial
end


function compute_random_utility(ξ_ij, mp = mp)

    u_ijt = mp.γ .- exp.(  -mp.ρ .* ξ_ij  ) 

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

function generate_price_transition_matrix(nStates)

    Π = rand(nStates,nStates)./100 # rand(J,nStates) Assume just one price state for now ...
    Π = Π./sum(Π, dims=2)

    return Π
end

function obtain_price_states_over_time(T, Π, inital_state = 3)

    # Obtain states (fake data) based on transition probabilities: 
    m_states = zeros(T)

    # Assume initial state = 0
    m_states[1] = inital_state

    for tt = 2:T
        m_states[tt] = rand(Categorical(Π[Int(m_states[Int(tt-1)]),:]))
    end

    return Int.(m_states)
end


function generate_model_simulation(mp, dims, cheb_params, quad_params, maxiter = 2000)

    # Unpack Parameters:
    α = mp.α
    β = mp.β
    σ_v = mp.σ_v
    σ_0 = mp.σ_0
    μ_0 = mp.μ_0
    ϑ = mp.ϑ

    # Unpack Dimensions:
    J = dims.J                  # Number of products
    nStates = dims.nStates      # Number of price states
    T = dims.T            # Time periods for updating rule

    # Generate fake prices, transition matrix and sequence of states:
    P =  [0.217756  0.0969838;
            0.774392  0.582034;
            0.662304  0.826249;
            0.690314  0.343262;
            0.150091  0.439394]

    Π = generate_price_transition_matrix(nStates)

    m_states = obtain_price_states_over_time(T, Π) 

    # Generate logit shocks:
    ϵ_jt = rand(Gumbel(0,1), J, T)

    # Get nodes and weights from quadrature assuming only state is quality:
    q_nodes, q_weights = GaussHermite(quad_params.N_dim, quad_params.N_States)

    # Generate arrays where we will store simulations: 
    
    # 1. Prior mean and variance Placeholders:
    μ_prior = zeros(J, T) ; μ_prior[:,1] .= μ_0
    σ_prior = zeros(J, T) ; σ_prior[:,1] .= σ_0

    # 2. Utility Placeholder:
    u_prime_sj = zeros(nStates, J)

    # 3. Choice:
    D_jt = zeros(J, T)
    
    V_all = zeros(J, T)

    # Add possibility to choose outside option.
    for tt = 1:(T-1)

        Vdiff = 10

        EV = randn(nStates,2)
        
        EV_next = zeros(nStates,2)    

        its = 1

        # COMPUTING VALUE FUNCTION: (add more states for prices, but more importantly change contraction ...)
        while (Vdiff > 0.00000001) & (its < maxiter)

            # Compute beliefs about tomorrow's quality experience:
            ξ_prime_ij = q_nodes .* σ_prior[:,tt]' .+ μ_prior[:,tt]'

            for jj = 1:J
                
                #2. Apply Chebyshev's approximation to tomorrow's random utility:
                # Initialize Chebyshev approximator
                cheb = initializeChebyshevApproximator(cheb_params.N_dim, 
                                                        cheb_params.N_degree, 
                                                        cheb_params.N_nodes , 
                                                        [minimum(ξ_prime_ij[:,jj])], 
                                                        [maximum(ξ_prime_ij[:,jj])]
                                                        )

                # Calculate Chebyshev regression coefficients to approximate f
                cheb = calculateChebyshevCoefficients(compute_random_utility, cheb)

                # Approximate random utility using cheby polynomials (First case, only) (n_nodes times n_price_states)
                
                u_prime_sj[:,jj] = evaluateChebyshev(ξ_prime_ij[:,jj], cheb)' * q_weights .- α .* P[:,jj]

                EV_next[:,jj] = (u_prime_sj[:,jj] .+ β .* EV[:,jj])' * Π' # TODO: Generalize this ...

                EV_next[:,jj] = (u_prime_sj[:,jj] .+ β .* EV[:,jj])' * Π'
    
            end

            Vdiff = maximum(abs.(EV_next.-EV))

            EV = copy(EV_next) # Important, copy, not input right away...

            print("Period number:",tt, "; V_diff = ", Vdiff,"\n")
            
            its += 1
        end

        # Draw quality shock at time t:
        v_ij = rand(Normal(0, σ_v.^2), (J,1))

        # Compute today's quality experience:
        ξ_jt = v_ij + ϑ'
        
        # Compute value function based on current shocks and EV beliefs:
        flow_utility = compute_random_utility(ξ_jt) .- α .* P[m_states[tt],:] .+ ϵ_jt[:,tt]

        V_all[:,tt] = flow_utility + β.* EV[m_states[tt],:]
        
        D_jt[:,tt] = V_all[:,tt] .== maximum(V_all[:,tt])

        # Update Beta:
        β_update = update_beta(σ_prior[:,tt], σ_v)

        # Update Prior Mean for next period: 
        μ_prior[:, tt+1] = update_quality_mean(μ_prior[:,tt], D_jt[:,tt], v_ij, β_update) 

        # Update Prior Variance for next period: 
        σ_prior[:, tt+1] = update_quality_variance(σ_prior[:,1], σ_prior[:,tt], D_jt[:,1:tt]) 

        T = tt-1
    end
    
    return V_all[:,1:T], D_jt[:,1:T], μ_prior[:,1:T], σ_prior[:,1:T]
end



dims = ModelDimensions(1000, # nObs
                        2,   # J
                        5,   # nStates
                        10,  # nPolynomials
                        20   # T
                        )  

mp   = ModelParameters(.9,      # ρ
                        .79,    # α
                        .5,     # γ
                        .995,   # β
                        .5,     # σ_v
                        .5,     # σ_0
                        1,      # μ_0
                        [.5 .5] # ϑ
                        )

cheb_params = ChebyshevParameters(  1, # N_dim
                                    5, # N_degree
                                    9  # N_nodes
                                    )

quad_params = QuadratureParameters( 1, # N_dim
                                    50, # N_states
                                    )

# Model 1: β = 0,998; 0,995; 0,99; 0;
V_all, D_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params)






































function generate_model_simulation(mp, dims, cheb_params, quad_params, maxiter = 200)

    # Unpack Parameters:
    α = mp.α
    β = mp.β
    σ_v = mp.σ_v
    σ_0 = mp.σ_0
    μ_0 = mp.μ_0
    ϑ = mp.ϑ

    # Unpack Dimensions:
    J = dims.J
    nStates = dims.nStates
    T = maxiter


    # Generate fake prices:
    P = rand(J,1) # rand(J,nStates) Assume just one price state for now ...

    # Generate logit shocks:
    ϵ_jt = rand(Gumbel(0,1), J, T)

    # Get nodes and weights from quadrature assuming only state is quality:
    q_nodes, q_weights = GaussHermite(quad_params.N_dim, quad_params.N_States)

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

    # COMPUTING VALUE FUNCTION: (add more states for prices, but more importantly change contraction ...)
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
            cheb = initializeChebyshevApproximator(cheb_params.N_dim, 
                                                    cheb_params.N_degree, 
                                                    cheb_params.N_nodes , 
                                                    [minimum(ξ_prime_ij[:,jj])], 
                                                    [maximum(ξ_prime_ij[:,jj])]
                                                    )

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

    T = tt-1
    
    return V_all[:,1:T], D_jt[:,1:T], μ_prior[:,1:T], σ_prior[:,1:T]
end
