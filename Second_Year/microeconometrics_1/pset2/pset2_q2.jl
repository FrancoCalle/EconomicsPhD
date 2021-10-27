using Pkg

#Load packages ...

using FixedEffectModels
using Distributions
using LinearAlgebra
# using StatsBase
# using CategoricalArrays
using DataFrames
using Plots
using StatFiles
using Chain
using Tables
# using CSV
using Optim
using Random
using StatsPlots

include(joinpath(@__DIR__,"..", "..", "..","fc_toolkit.jl"))

# Simple heteroskedastic case:
#Set sample size 
n=1000 
#Generate a single covariate. 
#I’m assuming the covariate follows a uniform[0, 100] distribution. 
z = rand(Normal(0,1),n)
#Generate an error term, which is normally distributed with mean 0 and SD 4 
sd = 0.8 #rand(Uniform(0,4),n)
v = vcat(rand.(Normal.(0, sd.* sqrt.(z.^2)),1)...)
#Choose parameters (intercept and slope) and generate the outcome: 
β = 1
x = β .* z + v

plot(z,x,seriestype = :scatter)




# Heteroskedastic IV case:
#Set sample size 

function heteroskedastic_iv(ρ_uv, γ, β, n)

    #Generate a single covariate. I assume it comes from a NoncentralChisq 
    Z = rand(NoncentralChisq(2,8), n)
    
    #Generate error terms U and V which are normally distributed with mean 0
    #Compute sigma
    Σ = Matrix{Float64}(I, 2, 2); Σ[1,2] = ρ_uv; Σ[2,1] = ρ_uv

    #Heteroskedastic errors
    vu_homoskedastic = hcat([rand(MvNormal(zeros(2), Σ * z), 1) for z in Z]...)'
    #Homoskedac errors
    vu_heteroskedastic = hcat([rand(MvNormal(zeros(2), Σ * 1), 1) for z in Z]...)'

    #Choose parameters (intercept and slope) and generate the outcome: 
    function dgp(vu)
        
        v = vu[:,1]
        u = vu[:,2]
        x = γ .* Z + v
        y = β .* x + u
                
        return y, x
    end

    y_hom, x_hom = dgp(vu_homoskedastic)
    y_het, x_het = dgp(vu_heteroskedastic)

    dgp_variables = Dict(y_hom => y_hom, 
            x_hom => x_hom,
            y_het => y_het,
            x_het => x_het)

    v = vu_heteroskedastic[:,1]
    u = vu_heteroskedastic[:,2]

    iv_shift = mean(Z.^2 .* u .* v)/mean(Z.^2 .* v.^2)
    ols_bias = mean(u .* v)/mean(v.^2)        

    return iv_shift, ols_bias, dgp_variables

end



function simulate(ρ_uv=0.1, γ = 0.02, β = 0, n=1000,M = 5000)

    iv_shift, ols_bias, dgp_variables = heteroskedastic_iv(ρ_uv, γ, β, n);

    
    bias_results = zeros(M,2)
    bias_distance = zeros(M,1)

    for m in 1:M
        iv_shift, ols_bias, dgp_variables = heteroskedastic_iv(ρ_uv, γ, β, n)
        bias_results[m,:] = [iv_shift, ols_bias]
        bias_distance[m] = iv_shift - ols_bias
    end

    return bias_results, bias_distance
end


ρ_uv_list = [0.99 0.4 0.2 0.1 0.05]
bias_dict = Dict()
for ρ in ρ_uv_list
    bias_results, bias_distance = simulate(ρ)
    bias_dict[Symbol("bias_distance", ρ)] = bias_distance./ρ
    bias_dict[Symbol("bias_results", ρ)] = bias_results
end


# density(bias_dict[Symbol("bias_distance", 0.99)], linewidth = 2, color = :red, label = "ρ = "*string(0.4), linestyle=:dash)
density(bias_dict[Symbol("bias_distance", 0.4)], linewidth = 2, color = :red, label = "ρ = "*string(0.4), linestyle=:dash)
density!(bias_dict[Symbol("bias_distance", 0.2)], linewidth = 2, color = :red, label = "ρ = "*string(0.2), linestyle=:dot)
density!(bias_dict[Symbol("bias_distance", 0.1)], linewidth = 2, color = :red, label = "ρ = "*string(0.1), linestyle=:dashdot)
density!(bias_dict[Symbol("bias_distance", 0.05)], linewidth = 2, color = :red, label = "ρ = "*string(0.05) )
vline!([0],color=:gray, linestyle=:dash, label="",linewidth=1)
savefig("heteroskedasticity_bias.pdf")

histogram(bias_results[:,1])
histogram!(bias_results[:,2])
ρ = 0.1
