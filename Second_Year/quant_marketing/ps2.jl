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


function compute_cost(v0,v1)
    
    return v1 .- v0 .- (v1[1] - v0[1])  

end

# Execute code:...


β = .98
δ = .2
r = 1
N = 10
nsim=1000

v0, v1 = simulation(β, δ, r, N, nsim)

c = (v1 .- v0) .- (v1[1] .- v0[1])



# Plotting figures: Changes in Type and Discount Rates:

delta_list = Array(0.01:0.01:1)
cost_placeholder = zeros(size(delta_list,1), N, 4)

jj=1

for delta in delta_list

    v0,v1 = simulation(.99, delta, 1, 10)
    cost_placeholder[jj,:,1] = compute_cost(v0,v1)

    v0,v1 = simulation(.97, delta, 1, 10)
    cost_placeholder[jj,:,2] = compute_cost(v0,v1)

    v0,v1 = simulation(.95, delta, 1, 10)
    cost_placeholder[jj,:,3] = compute_cost(v0,v1)

    v0,v1 = simulation(.8, delta, 1, 10)
    cost_placeholder[jj,:,4] = compute_cost(v0,v1)

    jj += 1

end

scatter(cost_placeholder[2:100, 7, 1], shape=:+, color=:gray, label="β: .99")
scatter!(cost_placeholder[2:100, 7, 2], shape=:o, color=:gray, label="β: .97")
scatter!(cost_placeholder[2:100, 7, 3], shape=:x, color=:gray, label="β: .95")
savefig("PS2_Plot_values_type_discount.pdf")



# Plotting figures: Changes reward utility:

delta_list = Array(0.01:0.01:1)
cost_placeholder_v2 = zeros(size(delta_list,1), N, 3)

jj=1

for delta in delta_list

    v0,v1 = simulation(.99, delta, 1, 10)
    cost_placeholder_v2[jj,:,1] = compute_cost(v0,v1)

    v0,v1 = simulation(.99, delta, 2, 10)
    cost_placeholder_v2[jj,:,2] = compute_cost(v0,v1)

    v0,v1 = simulation(.99, delta, 3, 10)
    cost_placeholder_v2[jj,:,3] = compute_cost(v0,v1)

    jj += 1

end


scatter(cost_placeholder_v2[2:100, 7, 1], shape=:+, color=:gray, label="r: 1")
scatter!(cost_placeholder_v2[2:100, 7, 2], shape=:o, color=:gray, label="r: 2")
scatter!(cost_placeholder_v2[2:100, 7, 3], shape=:x, color=:gray, label="r: 3")
savefig("PS2_plot_different_reward_utility.pdf")



# Part 5: 

# Get Median value of switch cost

delta_list = Array(0.05:0.05:5)
cost_placeholder = zeros(size(delta_list,1), N, 4)

jj=1

for delta in delta_list

    v0,v1 = simulation(.99, delta, 1, 10)
    cost_placeholder[jj,:,1] = compute_cost(v0,v1)

    v0,v1 = simulation(.97, delta, 1, 10)
    cost_placeholder[jj,:,2] = compute_cost(v0,v1)

    v0,v1 = simulation(.95, delta, 1, 10)
    cost_placeholder[jj,:,3] = compute_cost(v0,v1)

    v0,v1 = simulation(.8, delta, 1, 10)
    cost_placeholder[jj,:,4] = compute_cost(v0,v1)

    jj += 1

end

# Median:
c1 = median(cost_placeholder[:,:,1], dims=2)[2:100]
c2 = median(cost_placeholder[:,:,2], dims=2)[2:100]
c3 = median(cost_placeholder[:,:,3], dims=2)[2:100]
c4 = median(cost_placeholder[:,:,4], dims=2)[2:100]

scatter(c1, shape=:+, color=:gray, label="β: .99")
scatter(c2, shape=:+, color=:gray, label="β: .99")
scatter(c3, shape=:+, color=:gray, label="β: .99")
scatter(c4, shape=:+, color=:gray, label="β: .99")


# Mean:
c1 = mean(cost_placeholder[:,2:end,1], dims=2)[2:100]
c2 = mean(cost_placeholder[:,2:end,2], dims=2)[2:100]
c3 = mean(cost_placeholder[:,2:end,3], dims=2)[2:100]
c4 = mean(cost_placeholder[:,2:end,4], dims=2)[2:100]

scatter(c1, shape=:+, color=:gray, label="β: .99")
scatter(c2, shape=:+, color=:gray, label="β: .99")
scatter(c3, shape=:+, color=:gray, label="β: .99")
scatter(c4, shape=:+, color=:gray, label="β: .99")

# Get ratio of switching costs






####### FIX:

# temp <- t(rep(0,20))

# delta_vec <- seq.int(0,5, 0.05)

# type_list <- c()

# for (del in delta_vec){
#   temp_result <- experiment(0.99, del, 3, 20)
#   temp_v0 = temp_result$v0
#   temp_v1 = temp_result$v1
#   temp_c = temp_v1 - temp_v0 - (temp_v1[1] - temp_v0[1])
  
#   type_list <- c(type_list, mean(exp(temp_v1 - temp_v0)))
  
#   temp <- rbind(temp, t(temp_c))
# }


# plot(type_list)


