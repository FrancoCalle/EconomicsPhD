using Pkg

#Install ...
# Pkg.activate(".") 
# Pkg.add("Distributions")
# Pkg.add("StatsBase")
# Pkg.add(["DataFrames","DataFramesMeta","Chain"])
# Pkg.add("Plots")
# Pkg.add("CategoricalArrays")
# Pkg.add("StatFiles")
# Pkg.add("Tables")
# Pkg.add("CSV")
# Pkg.add("Optim")
# Pkg.add("StatsPlots")
# Pkg.add("FixedEffectModels")
# Pkg.add("BenchmarkTools")

Pkg.instantiate()

#Load packages ...

using FixedEffectModels
using Distributions
using LinearAlgebra
using DataFrames
using Plots
using Optim
using Random
using Tables
using StatFiles
using StatsBase
using CategoricalArrays
using CSV
using Base.Threads
using BenchmarkTools

include(joinpath(@__DIR__,"..", "..", "..","fc_toolkit.jl"))


# Start question 3:

struct define_parameters

    N::Int64
    T::Int64
    α::Float64
    β::Float64
    θ::Float64
    ρ::Float64

end

function data_generating_process(parameters,seed=MersenneTwister(1234))

    N = parameters.N
    T = parameters.T
    α = parameters.α
    β = parameters.β
    θ = parameters.θ
    ρ = parameters.ρ

    ID_i  = zeros(N, T+1)
    TT    = zeros(N, T+1)
    E_i   = zeros(N, T+1)
    ϵ_it  = zeros(N, T+1)
    U_it  = zeros(N, T+1)
    V_it  = zeros(N, T+1)

    treatment  = zeros(N, T+1)
    Y0_it = zeros(N, T+1)
    Y1_it = zeros(N, T+1)
    Y_it  = zeros(N, T+1)

    # ϵ_it[:,1] = rand(seed,Normal(),N)    
    E_i[:,1]  = sample(MersenneTwister(1234),2:5, N)

    for ii = 1:N
        
        for tt = 2:(T+1)

            TT[ii,tt] = tt-1

            ID_i[ii,tt] = ii

            E_i[ii,tt] = E_i[ii,tt-1] 

            treatment[ii,tt] =  TT[ii,tt]-E_i[ii,tt] >= 0

            ϵ_it[ii,tt] = rand(Normal())

            V_it[ii,tt] = rand(Normal())

            U_it[ii,tt] = ρ * U_it[ii,tt-1] + ϵ_it[ii,tt]

            Y0_it[ii,tt] = α + β * E_i[ii,tt] + U_it[ii,tt]

            Y1_it[ii,tt] = α + β * E_i[ii,tt] + sin(tt-1 - θ * E_i[ii,tt]) + U_it[ii,tt] + V_it[ii,tt]

            if treatment[ii, tt] == 1

                Y_it[ii,tt] = Y1_it[ii,tt]

            else

                Y_it[ii,tt] = Y0_it[ii,tt]
            
            end

        end
        
    end

    df = DataFrame( ID  = ID_i[:,2:T+1][:],
                    TT  = TT[:,2:T+1][:], 
                    EE  = E_i[:,2:T+1][:], 
                    ϵ   = ϵ_it[:,2:T+1][:], 
                    V   = V_it[:,2:T+1][:], 
                    U   = U_it[:,2:T+1][:], 
                    YY0 = Y0_it[:,2:T+1][:], 
                    YY1 = Y1_it[:,2:T+1][:],
                    YY  = Y_it[:,2:T+1][:],
                    Treatment  = treatment[:,2:T+1][:],
                    D   = (TT.-E_i)[:,2:T+1][:]
                    );

    return df
end


# Part B: Fixed effects estimation ...
function estimation(df, cl=false, vartype="hom")
    
    N = Integer(maximum(df[:,:ID]))
    T = Integer(maximum(df[:,:TT]))
    CohortFixedEffects          = Matrix(reduce(hcat, [df[:,:EE] .== tt for tt in unique(df[:,:EE])]))
    TimeFixedEffects            = Matrix(reduce(hcat, [df[:,:TT] .== tt for tt in unique(df[:,:TT])]))
    relativeTimeFixedEffects    = Matrix(reduce(hcat, [df[:,:D] .== tt for tt in [-3, -2, 0, 1, 2, 3]]))    
    AllFixedEffects = Matrix(reduce(hcat, [relativeTimeFixedEffects, 
                                                CohortFixedEffects[:,2:end], 
                                                TimeFixedEffects[:,2:end]]))
    Y = Array(df[:,:YY])

    if cl
        cluster = CategoricalArray(df[:, :ID])
        ols = olsRegression(Y, ones(N*T,1), AllFixedEffects, cluster)
    else
        ols = olsRegression(Y, ones(N*T,1), AllFixedEffects)
    end

    # Test under the null of true (just for D_1)
    results = inference(ols ,vartype, ones(size(ols.β)).*sin(1)) 

    return results.β, results.se, results.t, results.p, ols

end 



function monte_carlo(parameters, m = 12, α=0.025)
    
    parameters_placeholder = zeros(m,14)

    for mm in 1:m
        parameters_placeholder[mm,:],_ = estimation(data_generating_process(parameters,MersenneTwister())) 
    end    

    p_lower = [quantile(parameters_placeholder[:,jj], α)   for jj in 1:size(parameters_placeholder,2)]
    
    p_upper  = [quantile(parameters_placeholder[:,jj], 1-α) for jj in 1:size(parameters_placeholder,2)]
    
    avg = mean.(eachcol(parameters_placeholder))

    return avg, p_upper, p_lower

end


function average_treatment_effect(df)

    t1 = mean(df[(df.TT .== 2).& (df.EE .>= 3), :YY]) - mean(df[(df.TT .== 1).& (df.EE .>= 3), :YY])
    t2 = mean(df[(df.TT .== 3).& (df.EE .>= 4), :YY]) - mean(df[(df.TT .== 2).& (df.EE .>= 4), :YY])
    ATE2 = (mean(df[(df.TT .== 3) .& (df.EE .== 2), :YY]) 
                - mean(df[(df.TT .== 1) .& (df.EE .== 2), :YY]) 
                - t1 
                - t2)

    return ATE2
end


function monte_carlo_ate(parameters, m = 12, α=0.025)
    
    parameters_placeholder = zeros(m)

    for mm in 1:m
        parameters_placeholder[mm] = average_treatment_effect(data_generating_process(parameters,MersenneTwister())) 
    end    

    p_lower = quantile(parameters_placeholder, α)
    
    p_upper  = quantile(parameters_placeholder, 1-α) 
    
    avg = mean(parameters_placeholder)

    return avg, p_upper, p_lower

end


function monte_carlo_size_test(parameters, cl = false, vartype="hom", m = 12, α=0.05)
    
    parameters_placeholder = zeros(m,14)

    pval_placeholder = zeros(m,14)

    for mm in 1:m
        parameters_placeholder[mm,:],
        _,
        _,
        pval_placeholder[mm,:],
        _ = estimation(data_generating_process(parameters,MersenneTwister()), cl, vartype) 
    end    

    # Apply test for all draw results:
    power = mean(pval_placeholder[:,5] .<α)

    return power

end



# A and B: Simulation θ = -2... 
true_parameters = define_parameters(1000, 5, -0.2, 0.5, -2, 0.5)
avg1k, p_up1k, p_low1k = monte_carlo(true_parameters,  50, 0.05);

true_parameters = define_parameters(10000, 5, -0.2, 0.5, -2, 0.5)
avg10k, p_up10k, p_low10k = monte_carlo(true_parameters,  50, 0.05);

xAxis = [-3, -2, 0, 1, 2, 3]
plot((xAxis, avg10k[2:7]), markershape = :square, label=["N=10K, ci=0.95"], fillrange = [p_up10k[2:7]], fillalpha=0.3, c=:blue)
plot!((xAxis, avg10k[2:7]), markershape = :square, label=nothing, fillrange = [p_low10k[2:7]], fillalpha=0.3, c=:blue)
plot!((xAxis, avg1k[2:7]), markershape = :square, label=["N=1K, ci=0.95"], fillrange = [p_up1k[2:7]], fillalpha=0.3, c=:orange)
plot!((xAxis, avg1k[2:7]), markershape = :square, label=nothing, fillrange = [p_low1k[2:7]], fillalpha=0.3, c=:orange)
plot!(legend=:bottomleft)
savefig("Q3_Pb_Theta-2.pdf")


# C: Simulation θ = 0...
true_parameters = define_parameters(1000, 5, -0.2, 0.5, 0, 0.5)
avg1k, p_up1k, p_low1k = monte_carlo(true_parameters,  50, 0.05);

true_parameters = define_parameters(10000, 5, -0.2, 0.5, 0, 0.5)
avg10k, p_up10k, p_low10k = monte_carlo(true_parameters,  50, 0.05);

xAxis = [-3, -2, 0, 1, 2, 3]
plot((xAxis, avg10k[2:7]), markershape = :square, label=["N=10K, ci=0.95"], fillrange = [p_up10k[2:7]], fillalpha=0.3, c=:blue)
plot!((xAxis, avg10k[2:7]), markershape = :square, label=nothing, fillrange = [p_low10k[2:7]], fillalpha=0.3, c=:blue)
plot!((xAxis, avg1k[2:7]), markershape = :square, label=["N=1K, ci=0.95"], fillrange = [p_up1k[2:7]], fillalpha=0.3, c=:orange)
plot!((xAxis, avg1k[2:7]), markershape = :square, label=nothing, fillrange = [p_low1k[2:7]], fillalpha=0.3, c=:orange)
plot!(legend=:bottomleft)
savefig("Q3_Pc_Theta-0.pdf")


# C: Simulation θ = 1...
true_parameters = define_parameters(1000, 5, -0.2, 0.5, 1, 0.5)
avg1k, p_up1k, p_low1k = monte_carlo(true_parameters,  50, 0.05);

true_parameters = define_parameters(10000, 5, -0.2, 0.5, 1, 0.5)
avg10k, p_up10k, p_low10k = monte_carlo(true_parameters,  50, 0.05);

xAxis = [-3, -2, 0, 1, 2, 3]
plot((xAxis, avg10k[2:7]), markershape = :square, label=["N=10K, ci=0.95"], fillrange = [p_up10k[2:7]], fillalpha=0.3, c=:blue)
plot!((xAxis, avg10k[2:7]), markershape = :square, label=nothing, fillrange = [p_low10k[2:7]], fillalpha=0.3, c=:blue)
plot!((xAxis, avg1k[2:7]), markershape = :square, label=["N=1K, ci=0.95"], fillrange = [p_up1k[2:7]], fillalpha=0.3, c=:orange)
plot!((xAxis, avg1k[2:7]), markershape = :square, label=nothing, fillrange = [p_low1k[2:7]], fillalpha=0.3, c=:orange)
plot!(legend=:bottomleft)
savefig("Q3_Pc_Theta-1.pdf")


# Part D: Montecarlo for ATE
ate = zeros(3); ate_up = zeros(3); ate_lo = zeros(3); ate_true = zeros(3)
θ_list = [-2,0,1]

for ii in 1:3
    θ = θ_list[ii]
    true_parameters = define_parameters(1000, 5, -0.2, 0.5, θ, 0.5)
    ate[ii], ate_up[ii], ate_lo[ii] = monte_carlo_ate(true_parameters,  1000, 0.05);
    ate_true[ii] = sin(3 - 2 * θ)
end

scatter(ate_true, ate, markershape = :square, label="Coefficients", c=:blue, xlabel = "True Parameters", ylabel = "Simulation Average")
plot!([0,1],[0,1] , label="45 degree line")
scatter!(legend=:topleft)
savefig("ATE_estimation.pdf")



# Part E: Monte Carlo for to get estimator power...

# 1. Classical asymptotic variance estimator under homoskedacity

n_list = [20, 50, 200]

ρ_list = [0, .5, 1]

size_results_1 = zeros((length(n_list),length(ρ_list)))

size_results_2 = zeros((length(n_list),length(ρ_list)))

size_results_3 = zeros((length(n_list),length(ρ_list)))

size_results_4 = zeros((length(n_list),length(ρ_list)))

# size_results_5 = zeros((length(n_list),length(ρ_list)))


for ii in 1:length(n_list)

    for jj in 1:length(ρ_list)

        true_parameters = define_parameters(n_list[ii], 5, -0.2, 0.5, 1, ρ_list[jj])

        # Specification 1
        size_results_1[jj,ii] = monte_carlo_size_test(true_parameters, false, "hom", 2000);

        # Specification 2: HC(1), heteroskedasticity proof
        size_results_2[jj,ii] = monte_carlo_size_test(true_parameters, false, "het", 2000);

        # Specification 3: Cluster-robust asymptotic variance estimator,
        size_results_3[jj,ii] = monte_carlo_size_test(true_parameters, true, "clust", 2000);

        # Specification 4: Wild bootstrap
        size_results_4[jj,ii] = monte_carlo_size_test(true_parameters, true, "wild", 2000);

    end

end


#n is rows
#ρ is columns

results_table = DataFrame(  NObs            = [20,20,20, 50,50,50, 200,200,200],
                            Rho             = [0, .5, 1, 0, .5, 1, 0, .5, 1],
                            Homoskedastic   = size_results_1[:],
                            Heteroskedastic = size_results_2[:], 
                            Cluster         = size_results_3[:], 
                            WildCluster     = size_results_4[:], 
                            );

println(results_table)

CSV.write("Q3_PE_size_test.csv",  round.(results_table,sigdigits = 3), writeheader=false)





