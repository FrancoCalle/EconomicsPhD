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
# using Base.Threads



function get_decision(mileage)

    d_t = zeros(size(mileage))

    for ii = 1:size(mileage,1)-1 
        if ii > 1 
            d_t[ii] = mileage[ii+1] < mileage[ii]  
        end
    end

    return d_t

end


function classify(x; c=bin_edges, K=10)
    
    s = 0
    for kk = 1:K
        if c[kk] <= x < c[kk+1]
            s = kk
        end
    end

    return s

end

function generate_transition_matrix(mileage, d_t, c; T=5000, K = 10)

    P_0 = zeros(K,K)
    P_1 = zeros(K,K)

    for tt = 1:T-1

        m_next = mileage[tt+1]
        m = mileage[tt]
        d = d_t[tt]
        
        if d == false
            ii = classify(m; c, K)
            jj = classify(m_next; c, K)
            P_0[ii,jj] += 1
        else
            ii = classify(m; c, K)
            jj = classify(m_next; c, K)
            P_1[ii,jj] += 1
        end

    end

    # Add random shock (cheating...)
    P_0 .+= rand(K,K)./100 
    P_1 .+= rand(K,K)./100

    Π0 = P_0./sum(P_0,dims=2)
    Π1 = P_1./sum(P_1,dims=2)

    return Π0, Π1

end


function flow_utility(x , d ; θ_1=.3, θ_2=.7, θ_3 = .011)

    if d == 0
        u = - θ_1 * x - θ_2 * (x/100)^2
    else
        u = -θ_3 
    end
    
    return u
end







# Execute code:...

mileage_data = DataFrame(CSV.File("ps2_ex3.csv"));

mileage = mileage_data[:milage]

d_t = get_decision(mileage)

# 3. Discretize mileage in K

T = size(mileage,1)

K = 10

bin_edges = Array(LinRange(minimum(mileage),maximum(mileage),K+1)); bin_edges[end]= bin_edges[end] +1 

x = 100

# Transition probabilities:

Π0, Π1 = generate_transition_matrix(mileage, d_t, bin_edges; T, K)

# Compute contraction:

function obtain_continuation_value()

    β=0.98
    mileage_state = classify.(mileage; K = K)
    y_discrete = [mean(mileage[mileage_state.==ms]) for ms = 1:K]

    U = zeros(K, 2)
    U[:,1] = flow_utility.(y_discrete,0)
    U[:,2] = flow_utility.(y_discrete,1)

    EV = randn(K,2)
    EV_next = zeros(K,2)

    tol = 10

    while tol > 10^-15

        EV_next[:,1] = log.(exp.(U[:,1] .+ β .* EV[:,1]) .+ exp.(U[:,2] .+ β .* EV[:,2]))' * Π0'
        EV_next[:,2] = log.(exp.(U[:,1] .+ β .* EV[:,1]) .+ exp.(U[:,2] .+ β .* EV[:,2]))' * Π1'
        tol = maximum(abs.(EV_next.-EV))
        print(tol, "\n")
        EV = copy(EV_next) # Important, copy, not input right away...

    end

    return EV
end

# Check what is the regenerative property.

# 5. Derive the conditional choice probabilities using EV (x, d) and θ.
# 6. Reduce the state space of EV using the regenerative property.
# 7. Rewrite the fixed point equation as a matrix equation.
# 8. Write a function that solves the fixed-point equation using Rust’s algorithm.
# 9. Write a function that computes the likelihood of the sample for any parameter
# θ.
# 10. Estimate the model parameters θ. Use β = .999



