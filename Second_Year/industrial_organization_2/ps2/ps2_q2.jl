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
using RDatasets
using StatsPlots

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

        n_t = 1
        
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


function model_prediction(estimationData, mp; Eta_S, Epsilon_S, S=20)

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

function gmm_objective(parameters, I, N, estimationData; Eta_S, Epsilon_S, S)

    #Unpack data:
    T = estimationData.T;
    K = estimationData.K;
    X = estimationData.X; 
    Z = estimationData.Z;
    I = estimationData.I;
    N = estimationData.N';
    Market= estimationData.Market;


    mp0 = ModelParameters(parameters[1],
                            parameters[2],
                            parameters[3],
                            parameters[4],
                            1/(1+exp(parameters[5])))

    

    I_expected, N_expected = model_prediction(estimationData, mp0; Eta_S, Epsilon_S, S)

    ξ_micro = I .- I_expected

    ξ_macro = N .- N_expected

    obj_micro = Z[:]' * ξ_micro[:] * ξ_micro[:]' * Z[:]

    obj_macro = X[:]' * ξ_macro[:] * ξ_macro[:]' * X[:]

    obj_total = obj_micro + obj_macro 

    return obj_total

end


# Execute code:
#--------------

# Set model Parameters
mp = ModelParameters(1,2,6,3,0.8)

# Set model Data for DGP 
T = 10; K = 30;
X = exp.(rand(Normal(0,1), T));
Z = rand(Normal(0,1), K , T);
Firm = vcat([collect(1:30) for ii = 1:T]...);
Market = vcat([Int.(ones(K)*ii) for ii = 1:T]...);

md = ModelData(T, K, X, Z, Firm, Market);

# Compute DGP and obtain nash equilibrium
Π, I, N = generate_nash_equilibrium(md, mp)

# Gather all variables in data struct
estimationData = EstimationData(T, K, X, Z, I, N, Firm, Market)

# Estimation:
#-----------

# Get unobservable draws to feed the estimation procedure:
S = 100
Eta_S = [rand(Normal(0,1), T, S) for ss = 1:S];
Epsilon_S = [rand(Normal(0,1), T, K) for ss in 1:S];

# Grid Search:
#-------------

# Compute grid for α:
α_list = Array(-5:0.01:5)
obj_list_alpha = zeros(size(α_list,1))
for ii in 1:size(α_list,1)
    param_true = [α_list[ii],2,6,3,0.8]
    obj_list_alpha[ii] = gmm_objective(param_true, I, N, estimationData; Eta_S, Epsilon_S, S)
end

plot(α_list, obj_list_alpha, linewidth = 5, linecolor=:red, label="GMM Objective" )
vline!([1], linewidth=4, linecolor=:blue, label="α = 1")
savefig("q2_minimizing_at_alpha.pdf")


# Compute grid for β:
β_list = Array(-5:0.01:5)
obj_list_β = zeros(size(β_list,1))
for ii in 1:size(α_list,1)
    param_true = [1,β_list[ii],6,3,0.8]
    obj_list_β[ii] = gmm_objective(param_true, I, N, estimationData; Eta_S, Epsilon_S, S)
end


plot(β_list, obj_list_β, linewidth = 5, linecolor=:red, label="GMM Objective" )
vline!([2], linewidth=4, linecolor=:blue, label="β = 2")
savefig("q2_minimizing_at_beta.pdf")


# Compute grid for δ
δ_list = Array(4:0.01:9)
obj_list_δ = zeros(size(δ_list,1))
for ii in 1:size(δ_list,1)
    param_true = [1,2,δ_list[ii],3,0.8]
    obj_list_δ[ii] = gmm_objective(param_true, I, N, estimationData; Eta_S, Epsilon_S, S)
end

plot(δ_list, obj_list_δ, linewidth = 5, linecolor=:red, label="GMM Objective" )
vline!([6], linewidth=4, linecolor=:blue, label="δ = 6")
savefig("q2_minimizing_at_delta.pdf")

# Compute grid for γ:

γ_list = Array(1:0.01:5)
obj_list_γ = zeros(size(γ_list,1))
for ii in 1:size(γ_list,1)
    param_true = [1,2,6,γ_list[ii],0.8]
    obj_list_γ[ii] = gmm_objective(param_true, I, N, estimationData; Eta_S, Epsilon_S, S)
end

plot(γ_list, obj_list_γ, linewidth = 5, linecolor=:red, label="GMM Objective" )
vline!([3], linewidth=4, linecolor=:blue, label="γ = 3")
savefig("q2_minimizing_at_gamma.pdf")


# Compute grid for ρ:

ρ_list = Array(-10:0.01:10)
obj_list_ρ = zeros(size(ρ_list,1))
for ii in 1:size(ρ_list,1)
    param_true = [1,2,6,3,ρ_list[ii]]
    obj_list_ρ[ii] = gmm_objective(param_true, I, N, estimationData; Eta_S, Epsilon_S, S)
end

plot(1 ./(1 .+ exp.(ρ_list)), obj_list_ρ, linewidth = 5, linecolor=:red, label="GMM Objective" )
vline!([.8], linewidth=4, linecolor=:blue, label="ρ = .8")
savefig("q2_minimizing_at_rho.pdf")

# Report results for different initial conditions
#------------------------------------------------

nInit = 500
nParams = 5
shock = rand(nInit,nParams)
param_init_different = [1,2,6,3,log(1/0.8) - 1]' .+ shock .- mean(shock,dims=1)
param_init_results = zeros(nInit,nParams)
# Define Anonymous Function
func_anon(params) =  gmm_objective(params, I, N, estimationData; Eta_S, Epsilon_S, S)

for ii = 1:nInit
    # Proceed to parameter estimation:
    result = optimize(func_anon, param_init_different[ii,:], NelderMead(), 
                        Optim.Options(outer_iterations = 10000,
                                            iterations= 10000,
                                            show_trace=true,
                                            show_every=100
                                            )
                        )

    param_hat = Optim.minimizer(result)
    param_init_results[ii,:] = param_hat

end

param_init_results[:,end] = 1 ./(1 .+ exp.(param_init_results[:,end]))

param_init_different[:,end] = 1 ./(1 .+ exp.(param_init_different[:,end]))


# TODO: Add additional constraint to restrict values of \rho be less that 1 in absolute value...
using StatPlots
using PyPlot

pyplot()

# marginalhist(param_init_different[:,1], param_init_results[:,1], fc = :plasma)

# dataframe = DataFrame(initial_value = param_init_different[:,1], 
#                result = param_init_results[:,1]
#                )

# @df dataframe marginalhist(:initial_value, 
#                             :result, 
#                             nbins=30)






#--------------------------------------------------------


layout = @layout [a _
        b{0.8w,0.8h} c]

# Alpha:
default(fillcolor = :lightgrey, markercolor = :white, grid = false, legend = false)
plot(layout = layout, link = :both, size = (500, 500), margin = 5Plots.px)
scatter!(param_init_different[:,1], 
            param_init_results[:,1], 
            subplot = 2, 
            framestyle = :box            )
histogram!([param_init_different[:,1] param_init_results[:,1]], 
            subplot = [1 3], 
            orientation = [:v :h], 
            framestyle = :none,
            )

savefig("q2_2histogram_alpha.pdf")

# Beta:
default(fillcolor = :lightgrey, markercolor = :white, grid = false, legend = false)
plot(layout = layout, link = :both, size = (500, 500), margin = 5Plots.px)
scatter!(param_init_different[:,2], 
            param_init_results[:,2], 
            subplot = 2, 
            framestyle = :box            )
histogram!([param_init_different[:,2] param_init_results[:,2]], 
            subplot = [1 3], 
            orientation = [:v :h], 
            framestyle = :none,
            )

savefig("q2_2histogram_beta.pdf")


# Delta:
default(fillcolor = :lightgrey, markercolor = :white, grid = false, legend = false)
plot(layout = layout, link = :both, size = (500, 500), margin = 5Plots.px)
scatter!(param_init_different[:,3], 
            param_init_results[:,3], 
            subplot = 2, 
            framestyle = :box            )
histogram!([param_init_different[:,3] param_init_results[:,3]], 
            subplot = [1 3], 
            orientation = [:v :h], 
            framestyle = :none,
            )

savefig("q2_2histogram_delta.pdf")


# Gamma:
default(fillcolor = :lightgrey, markercolor = :white, grid = false, legend = false)
plot(layout = layout, link = :both, size = (500, 500), margin = 5Plots.px)
scatter!(param_init_different[:,4], 
            param_init_results[:,4], 
            subplot = 2, 
            framestyle = :box            )
histogram!([param_init_different[:,4] param_init_results[:,4]], 
            subplot = [1 3], 
            orientation = [:v :h], 
            framestyle = :none,
            )

savefig("q2_2histogram_gamma.pdf")



# Rho:
default(fillcolor = :lightgrey, markercolor = :white, grid = false, legend = false)
plot(layout = layout, link = :both, size = (500, 500), margin = 5Plots.px)
scatter!(param_init_different[:,5], 
            param_init_results[:,5], 
            subplot = 2, 
            framestyle = :box            )
histogram!([param_init_different[:,5] param_init_results[:,5]], 
            subplot = [1 3], 
            orientation = [:v :h], 
            framestyle = :none,
            )

savefig("q2_2histogram_rho.pdf")
