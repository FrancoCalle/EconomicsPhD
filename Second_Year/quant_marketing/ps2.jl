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


function generate_transition_matrix(N = 10)

    Π0 = Matrix(I, N,N)
    Π1 = zeros(N,N)
    for ii = 2:N
        Π1[ii,ii-1] = 1
    end
    Π1[1,end] = 1

    return Int.(Π0), Int.(Π1)

end


function utility(a, x , ϵ; δ=.3, r=.7, β = 0.98, N = 10)

    redeem = zeros(N)

    for ii = 1:N
        if x[ii] == N
            redeem[ii] = 1
        end
    end
    
    u =  mean(a .* (δ .+ ϵ .+ redeem' .* r) .+ (1 .- a) .* ϵ, dims=1)

    return u
end


function fixed_point(P, EV_0, EV_1, u_0, u_1, β)

    fp = P * log.(exp.(u_0' .+ β .* EV_0) .+ exp.(u_1' .+ β .* EV_1))
    
    return fp
end


function fixed_point_iteration(x0, 
                                β, 
                                δ, 
                                r, 
                                N, 
                                Π0, 
                                Π1, 
                                ϵ; 
                                EV_init = zeros(N), 
                                tol = 1e-10, 
                                max_iter=10000)

    EV_0 = copy(EV_init)
    EV_1 = copy(EV_init)

    u_0 = utility(0, x0, ϵ ; δ, r, β, N = 10)
    u_1 = utility(1, x0, ϵ ; δ, r, β, N = 10)

    for ii = 1:max_iter

        EV_0_prime = copy(EV_0)
        EV_1_prime = copy(EV_1)

        EV_0 = fixed_point(inv(Π0), EV_0, EV_1, u_0, u_1, β)
        EV_1 = fixed_point(inv(Π1), EV_0, EV_1, u_0, u_1, β)

        if sum(abs.(EV_0_prime .- EV_0 .+ EV_1_prime .- EV_1)) < tol
            print("Number of iterations: ", ii, "\n")
            break
        end

    end

    return EV_1, EV_0

end


function simulation(β, δ, r, N, nsim=1000)

    ϵ = rand(Gumbel(), nsim)
    x0 = 1:N
    Π0, Π1 = generate_transition_matrix(N)
    EV_1, EV_0 = fixed_point_iteration(x0,
                                        β, 
                                        δ, 
                                        r, 
                                        N, 
                                        Π0, 
                                        Π1, 
                                        ϵ)

    v0 = utility(0, x0, ϵ; δ, r, β, N)' .+ β * EV_0
    v1 = utility(1, x0, ϵ; δ, r, β, N)' .+ β * EV_1
    
    return v0, v1
end



# Execute code:...


β = .98
δ = .2
r = 1
N = 10
nsim=1000

v0, v1 = simulation(β, δ, r, N, nsim)

c = (v1 .- v0) .- (v1[1] .- v0[1])



# Plotting figures:

delta_list = Array(0.01:0.01:1)
temp = zeros(size(delta_list,1),N)

jj=1
for delta in delta_list

    v0,v1 = simulation(β, delta, 1, 10)
    cost = v1 .- v0 .- (v1[1] - v0[1])  
    temp[jj,:] = cost'
    jj += 1

end

plot(temp[2:100, 7])

