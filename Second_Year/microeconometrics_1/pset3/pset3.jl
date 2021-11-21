using Pkg

#Install ...
Pkg.activate(".") 
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add(["DataFrames","DataFramesMeta","Chain"])
Pkg.add("Plots")
Pkg.add("CategoricalArrays")
Pkg.add("StatFiles")
Pkg.add("Tables")
Pkg.add("CSV")
Pkg.add("Optim")
Pkg.add("StatsPlots")
Pkg.add("FixedEffectModels")
Pkg.add("BenchmarkTools")

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

    ID_i  = Matrix{Float64}(undef, N, T+1)
    TT    = Matrix{Float64}(undef, N, T+1)
    E_i   = Matrix{Float64}(undef, N, T+1)
    ϵ_it  = Matrix{Float64}(undef, N, T+1)
    U_it  = Matrix{Float64}(undef, N, T+1)
    V_it  = Matrix{Float64}(undef, N, T+1)

    treatment  = Matrix{Float64}(undef, N, T+1)
    Y0_it = Matrix{Float64}(undef, N, T+1)
    Y1_it = Matrix{Float64}(undef, N, T+1)
    Y_it  = Matrix{Float64}(undef, N, T+1)

    E_i[:,1] = round.(rand(MersenneTwister(1234),Uniform(.5, 5.5),N))

    for ii = 1:N
        
        for tt = 2:(T+1)

            TT[ii,tt] = tt-1

            ID_i[ii,tt] = ii

            E_i[ii,tt] = E_i[ii,tt-1] 

            treatment[ii,tt] =  TT[ii,tt]-E_i[ii,tt] >= 0

            ϵ_it[ii,tt] = rand(seed,Normal())

            V_it[ii,tt] = rand(seed,Normal())

            U_it[ii,tt] = ρ * U_it[ii,tt-1] + ϵ_it[ii,tt]

            Y0_it[ii,tt] = α + β * E_i[ii,tt] + U_it[ii,tt]

            Y1_it[ii,tt] = α + β * E_i[ii,tt] + sin(tt - θ * E_i[ii,tt]) + U_it[ii,tt] + V_it[ii,tt]

            if treatment[ii, tt] == 1

                Y_it[ii,tt] = Y1_it[ii,tt]

            else

                Y_it[ii,tt] = Y0_it[ii,tt]
            
            end

        end
        
    end

    df = DataFrame( ID  =ID_i[:,2:T+1][:],
                    TT  =TT[:,2:T+1][:], 
                    EE  =E_i[:,2:T+1][:], 
                    ϵ   =ϵ_it[:,2:T+1][:], 
                    V   =V_it[:,2:T+1][:], 
                    U   =U_it[:,2:T+1][:], 
                    YY0 =Y0_it[:,2:T+1][:], 
                    YY1 =Y1_it[:,2:T+1][:],
                    YY  =Y_it[:,2:T+1][:],
                    Treatment  = treatment[:,2:T+1][:],
                    D   = (TT.-E_i)[:,2:T+1][:]
                    );

    return df
end


# Part B: Fixed effects estimation ...
function estimation(df)
    
    N = Integer(maximum(df[:,:ID]))
    T = Integer(maximum(df[:,:TT]))
    depvar = [:YY]
    CohortFixedEffects          = Matrix(reduce(hcat, [df[:,:EE] .== tt for tt in unique(df[:,:EE])]))
    TimeFixedEffects            = Matrix(reduce(hcat, [df[:,:TT] .== tt for tt in unique(df[:,:TT])]))
    relativeTimeFixedEffects    = Matrix(reduce(hcat, [df[:,:D] .== tt for tt in [-3, -2, 0, 1, 2, 3, 4]]))    
    AllFixedEffects = Matrix(reduce(hcat, [relativeTimeFixedEffects, 
                                                CohortFixedEffects[:,1:end-1], 
                                                TimeFixedEffects[:,1:end-1]]))

    Y, _, _, _  = select_variables(df, depvar, [:D] )
    ols1 = olsRegression(Y, ones(N*T,1), AllFixedEffects)
    relative_time_fe = ols1.β[2:size(relativeTimeFixedEffects,2)+1]
    
    return relative_time_fe

end 


function monte_carlo(parameters, m = 12, α=0.025)
    
    parameters_placeholder = zeros(m,7)

    for mm in 1:m
        parameters_placeholder[mm,:] = estimation(data_generating_process(parameters)) 
    end    

    p_low = [quantile(parameters_placeholder[:,jj], α) for jj = 1:size(parameters_placeholder,2)]
    p_up = [quantile(parameters_placeholder[:,jj], 1-α) for jj = 1:size(parameters_placeholder,2)]
    avg = mean.(eachcol(parameters_placeholder))

    return avg, p_up, p_low

end




function parallel_monte_carlo(parameters, m = 12, α=0.025)
    
    parameters_placeholder = zeros(m,7)
    rnglist = [MersenneTwister() for i in 1:nthreads()]

    @threads for mm in 1:m
        parameters_placeholder[mm,:] = estimation(data_generating_process(parameters, rnglist[threadid()])) 
    end    

    p_low = [quantile(parameters_placeholder[:,jj], α)   for jj = 1:size(parameters_placeholder,2)]
    p_up  = [quantile(parameters_placeholder[:,jj], 1-α) for jj = 1:size(parameters_placeholder,2)]
    avg = mean.(eachcol(parameters_placeholder))

    return avg, p_up, p_low

end




true_parameters = define_parameters(10000, 5, -0.2, 0.5, -2, 0.5)
df = data_generating_process(true_parameters)

println("Number of threads active: " , string(nthreads()))



avg, p_up, p_low = parallel_monte_carlo(true_parameters,  10000, 0.025)



@btime avg, p_up, p_low = monte_carlo(true_parameters,  2000, 0.05);
@btime avg, p_up, p_low = parallel_monte_carlo(true_parameters,  1000, 0.05);


avg, p_up, p_low = parallel_monte_carlo(true_parameters,  500, 0.05);


