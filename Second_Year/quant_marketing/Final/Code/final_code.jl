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
    # P =  [  0.317756  0.0969838;
    #         0.474392  0.102034;
    #         0.512304  0.1526249;
    #         0.520314  0.203262;
    #         0.550091  0.2540
    #         ]

    P =  [  0.317756  0.0969838;
    0.317756  0.0969838;
    0.317756  0.0969838;
    0.317756  0.0969838;
    0.317756  0.0969838;
            ]
            
    # P = abs.(randn(nStates, J))

    Π = generate_price_transition_matrix(nStates)

    m_states = obtain_price_states_over_time(T, Π) 

    # Generate logit shocks:
    # ϵ_jt = rand(Gumbel(0,1), J, T)

    # Get nodes and weights from quadrature assuming only state is quality:
    q_nodes, q_weights = GaussHermite(quad_params.N_dim, quad_params.N_States)

    # Generate arrays where we will store simulations: 
    
    # 1. Prior mean and variance Placeholders:
    μ_prior = zeros(J, T) ; μ_prior[:,1] .= μ_0
    σ_prior = zeros(J, T) ; σ_prior[:,1] .= σ_0

    # 3. Choice:
    D_jt = zeros(J + 1, T)

    Pr_jt = zeros(J + 1, T)
    
    V_all = zeros(J, T)

    TT = 0

    function get_beliefs_for_all_states(μ_prior, σ_prior, σ_v, v_ij, q_nodes = q_nodes, dims = dims)

        """
        ξ_prime_ij is a Nodes × JJ + 1 (previous choice state) × JJ (option beliefs) arrays
        that contains the nodes for each option conditional on choosing jj or the outside option in the 
        previous iteration.
        """

        # Update Beta:
        β_update = update_beta(σ_prior, σ_v)

        ξ_prime_ij = zeros(length(q_nodes),dims.J+1, J)

        # # Get shocks for outside option:
        ξ_prime_ij[:, dims.J+1, :] = q_nodes .* σ_prior' .+ μ_prior'

        for jj = 1:dims.J
            
            # Define the choice for each state:
            D = zeros(dims.J)
            D[jj] = 1
            
            # Update Prior Mean for next period ( J + outside option) beliefs:
            μ_prior_prime = update_quality_mean(μ_prior, D, v_ij, β_update)

            # # Update Prior Variance for next period: 
            σ_prior_prime = update_quality_variance(σ_prior, σ_prior, D) 

            # # Get shocks for all periods
            ξ_prime_ij[:, jj, :] = q_nodes .* σ_prior_prime' .+ μ_prior_prime'
        
        end

        return ξ_prime_ij 
    end


    function get_expected_value(ξ_prime_all, EV, α = α, P = P, q_weights = q_weights, Π = Π)

        J = size(ξ_prime_all,3)
        
        expected_value = zeros(size(P,1), J+1, J)  ## First two entries are states, last one is choice:

        for jj = 1:J # Sum across all choices 

            for ss = 1:J+1

                # print(jj,ss, '\n')

                cheb = initializeChebyshevApproximator(cheb_params.N_dim, 
                                                        cheb_params.N_degree, 
                                                        cheb_params.N_nodes , 
                                                        [minimum(ξ_prime_all[:,ss,jj])], 
                                                        [maximum(ξ_prime_all[:,ss,jj])]
                                                        )

                # Calculate Chebyshev regression coefficients to approximate f
                cheb = calculateChebyshevCoefficients(compute_random_utility, cheb)

                util = evaluateChebyshev(ξ_prime_all[:,ss,jj], cheb)' * q_weights .- α .* P[:,jj]

                expected_value[:,ss,jj] = util .+ β .* Π *  log.(sum(exp.(EV[:,ss,:]), dims=2))

                
            end

        end     

        return expected_value
    end
    
    # TODO: Add possibility to choose outside option ...
    for tt = 1:(T-1)

        # 1. Draw Shocks: 

        # Draw quality shock at time t:
        v_ij = rand(Normal(0, σ_v.^2), (J,1))

        # Compute today's quality experience:
        ξ_jt = v_ij + ϑ'
    
        # Compute beliefs about tomorrow's quality experience (j , -j , outside opt):
        ξ_prime_all = get_beliefs_for_all_states(μ_prior[:,tt], σ_prior[:,tt], σ_v, v_ij)

        # Parameters for Iterations:
        Vdiff = 10

        EV = randn(size(P,1), J+1, J)
        
        EV_next = zeros(size(P,1), J+1, J) 

        its = 1


        # COMPUTING VALUE FUNCTION: 
        while (Vdiff > 0.00000001) & (its < maxiter)
                
            # Approximate random utility using cheby polynomials (First case, only) (n_nodes times n_price_states)
            EV_next = get_expected_value(ξ_prime_all, EV)   # Transition matrix same because of iid. (Otherwise it would be q_weights × Π )

            Vdiff = maximum(abs.(EV_next.-EV))

            EV = copy(EV_next) # Important, copy, not input right away...

            if its % 250 == 0.0
                print("Period number: ",tt, "; Iteration: ",its, "; V_diff = ", Vdiff,"\n")
            end
            
            its += 1
        end
        
        # Compute value function based on current shocks and EV beliefs:
        flow_utility = compute_random_utility(ξ_jt) .- α .* P[m_states[tt],:] #.+ ϵ_jt[:,tt]

        if tt == 1
            s_state = J+1
        else
            s_state = argmax(D_jt[:,tt-1])
        end

        V_all[:,tt] = flow_utility + β.* EV[m_states[tt], s_state,:]
        
        V_placeholder = V_all[:,tt] .- minimum(V_all[:,tt])

        Pr_jt[1:J,tt] = exp.(V_placeholder) ./ (1 .+ sum(exp.(V_placeholder),dims=1))

        Pr_jt[J+1,tt] = 1 - sum(Pr_jt[1:J,tt])

        D_jt[rand(Categorical(Pr_jt[:,tt])),tt] = 1

        # Update Beta:
        β_update = update_beta(σ_prior[:,tt], σ_v)

        # Update Prior Mean for next period: 
        μ_prior[:, tt+1] = update_quality_mean(μ_prior[:,tt], D_jt[1:J,tt], v_ij, β_update) 

        # Update Prior Variance for next period: 
        σ_prior[:, tt+1] = update_quality_variance(σ_prior[:,1], σ_prior[:,tt], D_jt[1:J,1:tt]) 

        TT = tt

    end
    
    return V_all[:,1:TT-1], D_jt[:,1:TT-1], Pr_jt[1:J,1:TT-1], μ_prior[:,1:TT-1], σ_prior[:,1:TT-1]

end



dims = ModelDimensions(1000, # nObs
                        2,   # J
                        5,   # nStates
                        10,  # nPolynomials
                        20   # T
                        )  

mp   = ModelParameters(.9,      # ρ
                        .5,    # α
                        .42,     # γ
                        .995,   # β
                        .5,     # σ_v
                        .5,     # σ_0
                        1,      # μ_0
                        [.5 .5 ] # ϑ
                        )


cheb_params = ChebyshevParameters(  1, # N_dim
                                    5, # N_degree
                                    9  # N_nodes
                                    )

quad_params = QuadratureParameters( 1, # N_dim
                                    50, # N_states
                                    )

# Model 1: β = 0,998; 0,995; 0,99; 0;
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)


# 3. Analyse and simulate the model

plot(V_all[1,:], label="Product 1")
plot!(V_all[2,:], label="Product 2")
savefig("Expected_value_function_over_time.pdf")

plot(Pr_jt[1,:], label="Product 1")
plot!(Pr_jt[2,:], label="Product 2")
savefig("ccp_over_time.pdf")

plot(Pr_jt[1,:], label="Product 1")
plot!(Pr_jt[2,:], label="Product 2")
savefig("ccp_over_time.pdf")

plot(μ_prior[1,:], label="Product 1")
plot!(μ_prior[2,:], label="Product 2")
plot!([0; 18], [.5; .5], lw=2, lc=:black, label="True Mean Quality")
savefig("average_quality_learning_parameter.pdf")

# Part b:
mean(D_jt, dims=2)

pr_jt = cumsum(D_jt, dims=2)./ cumsum(ones(T-2))' 

plot(pr_jt[1,:], label="Product 1")
plot!(pr_jt[2,:], label="Product 2")
savefig("probability_evolution_over_time.pdf")



# Question D:


dims = ModelDimensions(1000, # nObs
                        2,   # J
                        5,   # nStates
                        10,  # nPolynomials
                        30   # T
                        )  

mp   = ModelParameters(.9,      # ρ
                        .8,    # α
                        .48,     # γ
                        .995,   # β
                        .5,     # σ_v
                        .5,     # σ_0
                        1,      # μ_0
                        [.5 .5 ] # ϑ
                        )

V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)

# Question C: Probability evolution:
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 

plot(pr_jt[1,:], label="Product 1")
plot!(pr_jt[2,:], label="Product 2")
savefig("probability_evolution_over_time.pdf")

# Question D: Product switches (These eliminate the outside option, now we have to check how many switches there are):
D_jt_no_oo = D_jt[1:2,D_jt[3,:].!=1]
switches = zeros(1, size(D_jt_no_oo, 2))

for ii = 1:size(D_jt_no_oo, 2)

    if ii >1
        switches[1,ii] = argmax(D_jt_no_oo[:,ii-1]) == argmax(D_jt_no_oo[:,ii])
    end

end

switches_cumsum = cumsum(switches, dims=2)
switches_all = zeros(1, size(D_jt, 2))
switches_all[1,D_jt[3,:].!=1] .= switches_cumsum'

for ii = 1:size(switches_all,2)
    if ii >1 
        if switches_all[ii]==0
            switches_all[ii] = switches_all[ii-1]
        end
    end
end

inside_option_cumsum = cumsum(sum(D_jt[1:2,:], dims=1), dims=2)

plot((switches_all./inside_option_cumsum)'[3:end,1])
savefig("switches_ratio.pdf")


# Subsection C: Changes in probability evolution conditional on different prior means:

mp   = ModelParameters(.9,      # ρ
                        .8,     # α
                        .48,     # γ
                        .995,   # β
                        .5,     # σ_v
                        .5,     # σ_0
                        1,      # μ_0
                        [.5 .5 ] # ϑ
                        )


mp.μ_0 = .5
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot(pr_jt[1,:], label="μ = 0.5; ϑ = 0.5")

mp.μ_0 = 1
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot!(pr_jt[1,:], label=label="μ = 1; ϑ = 0.5")

mp.μ_0 = 2
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot!(pr_jt[1,:], label=label="μ = 2; ϑ = 0.5")

# Question C: Probability evolution:
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot(pr_jt[1,:], label=label="μ = 0.5; ϑ = 0.5")

savefig("probability_evolution_conditional_on_beliefs.pdf")

mp.μ_0 = .5

# Subsection E: Compare the simulations for different discount factors, 0:998; 0:995; 0:99; 0; and 
# for different levels of the risk aversion parameter : ... 

# Discount Parameters:
mp.β = .998

V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot(pr_jt[1,:], label="β = 0.998")

mp.β = .995
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot!(pr_jt[1,:], label="β = 0.995")

mp.β = .90
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot!(pr_jt[1,:], label="β = 0.90")

mp.β = .80
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot!(pr_jt[1,:], label="β = 0.80")

mp.β = .50
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot!(pr_jt[1,:], label="β = 0.80")

xlims!((5,28))
ylims!((0,.75))

savefig("different_discount_parameters.pdf")




# Risk Aversion Parameters:

mp.ρ = .9
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot(pr_jt[1,:], label="ρ = .9")

mp.ρ = .8
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot!(pr_jt[1,:], label="ρ = .8")

mp.ρ = .7
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot!(pr_jt[1,:], label="ρ = .7")

mp.ρ = .5
V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 
plot!(pr_jt[1,:], label="ρ = .5")

savefig("different_risk_aversion_parameters.pdf")


# Part A:
mp.β = .998
mp.ρ = .9

V_all, D_jt, Pr_jt, μ_prior, σ_prior = generate_model_simulation(mp, dims, cheb_params, quad_params, 5000)
pr_jt = cumsum(D_jt, dims=2) ./ cumsum(ones(size(D_jt,2)))' 

plot(V_all[1,:] .+ 2)
plot!(V_all[2,:] .+ 2)
savefig("value_function_simulation.pdf")
