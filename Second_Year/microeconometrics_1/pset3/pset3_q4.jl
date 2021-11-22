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

function load_dataset()

    cigs = DataFrame(CSV.File("adh-cigs.csv"))

    return cigs 

end


function diff_diff(df, controlState="IL")

    treatmentState = "CA"
    treatmentYear = 1988
    df.T = df.year .> treatmentYear

    treatedGroup = df[df.state .== treatmentState,:]
    controlGroup = df[df.state .== controlState,:]
    
    diffTreated = mean(treatedGroup[treatedGroup.T.==1,:cigs])-mean(treatedGroup[treatedGroup.T.==0,:cigs])  
    diffControl = mean(controlGroup[controlGroup.T.==1,:cigs])-mean(controlGroup[controlGroup.T.==0,:cigs])  
    
    ate = diffTreated - diffControl

    return ate
end

# Part 1: Get diff diff for each control state:
cigs = load_dataset()
controlStates = unique(cigs.state[cigs.state.!= "CA"])

ate_ss = zeros(length(controlStates))
for ii in 1:length(controlStates)
    ss = controlStates[ii]
    ate_ss[ii] = diff_diff(cigs,ss)
end

bar(controlStates,ate_ss,orientation = :horizontal, label=["Diff & Diff by State"])
savefig("diff_diff_by_state.pdf")


# Part 2: Compute synthetic control:

maxy = maximum(cigs.year); miny = minimum(cigs.year)

XControl = zeros(maxy-miny+1, length(controlStates))
for ii in 1:length(controlStates)
    ss = controlStates[ii]
    XControl[:,ii] = cigs[cigs.state.==ss,:cigs]
end
XTreated = cigs[cigs.state.=="CA",:cigs] ;


# Pre Policy Data - California
X1 = cigs[(cigs.state.=="CA") .& (cigs.year .<= 1988),:cigs] ;

# Pre Policy Data - Other States
X0 = XControl[1:(1988-miny+1),:]

# Define objective function to get synthetic data:
func(λ) = sqrt((X1-X0*λ)'*(X1-X0*λ));
λ₀ = ones(length(controlStates)).*.2


# Optimized objective function
res = optimize(func, λ₀, LBFGS(),Optim.Options(iterations = 1000))
λ_hat = Optim.minimizer(res)

# Use lambda to produce synthetic data for California:
plot(miny:maxy,XControl*λ_hat, label="Synthetic", linewidth=3, c=:gray)
plot!(miny:maxy,XTreated, label="True", linewidth=3, c=:red)
vline!([1988], linewidth=2, c=:green, linestyle=:dot, label=nothing) 
savefig("Synthetic_extrapolation_cigarretes.pdf")


# Compute Diff & Diff:
syntheticCal = XControl*λ_hat
ate = mean(XTreated[(1988-miny+1):end,1]) - mean(syntheticCal[(1988-miny+1):end])



