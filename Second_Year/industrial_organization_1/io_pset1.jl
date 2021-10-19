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
Pkg.instantiate()

#Load packages ...

using Distributions
using LinearAlgebra
using DataFrames
using StatFiles
# using Chain
using Tables
using CSV
using Optim
using Random


# Question 2

function production_function(k, l, α)
    
    # Define parameters:
    ρ = 0.8
    f_γ = Normal(0, 0.5)
    f_ϵ = Normal(0, 0.1)

    γ_i = rand(f_γ,1)[1]
    ϵ_it = rand(f_ϵ,1)[1]


    ω_it = ρ * ω_ilag + ϵ_it
    a_it = γ_i + ω_it
    y = a_it + α * l_it + (1-α) k_it

    return y
end


