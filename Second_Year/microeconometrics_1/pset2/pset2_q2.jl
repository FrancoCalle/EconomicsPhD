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
Pkg.instantiate()
Pkg.add("FixedEffectModels")

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





















