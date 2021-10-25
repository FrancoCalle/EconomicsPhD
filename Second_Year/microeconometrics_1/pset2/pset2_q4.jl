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
Pkg.add("Missings")
Pkg.add("FixedEffectModels")
Pkg.add("RCall")

#Load packages ...

using FixedEffectModels
using Distributions
using LinearAlgebra
using DataFrames
using Plots
using StatFiles
using Chain
using Tables # using CSV
using Optim
using Random
using StatsPlots
using CategoricalArrays
using Missings
using RCall 

include(joinpath(@__DIR__, "..", "..", "..","fc_toolkit.jl"))

census_data = DataFrame(load("matched_censusdata.dta"));

function setVariableNames()
  
  controlVariableList = ["base_hhpovrate0","prop_head_f_a0", "sexratio0_a", 
                            "prop_indianwhite0", "kms_to_road0", "kms_to_town0",  
                            "prop_matric_m0","prop_matric_f0", "baseline_hhdens0", 
                            "kms_to_subs0"]

  controlVariableNameList = ["Poverty Rate", "Female-headed HHs", "Adult sex ratio (F/M)", 
                                "Indian, white adults x 10", "Kilometers to road", "Kilometers to town", 
                                "Men with high school", "Women with high school", "Household density", 
                                "Kilometers from grid"]

  additionalControls = ["d_prop_waterclose", "d_prop_flush"]

  controlVariables = Dict(:controlVariableList => [Symbol(x) for x in controlVariableList], 
                          :controlVariableNameList => [Symbol(x) for x in controlVariableNameList],  
                          :additionalControls => [Symbol(x) for x in additionalControls])

  return controlVariables

end

controlVariables  = setVariableNames()


"""
'Function that generates table 3 (First Stage):
'
'@df is a dataframe which contains the variables used in regressions for table3
'@controlVariables is a list which contains lists of names and tag of variables  

set control variables:
"""->
function generateTable3(df = census_data, controlVariables)

  cluster = CategoricalArray(df[:, :placecode0])
  FixedEffects = Matrix(reduce(hcat, [df[:,:dccode0].==fe for fe in unique(df[:,:dccode0])]))

  #Specification 1:
  Y, X, _  = select_variables(df, [:T], [:mean_grad_new])
  fit = olsRegression(X ./10 , Y, nothing, cluster)
  inference(fit)

  #Specification 2:
  x_names = vcat([:mean_grad_new],controlVariables[:controlVariableList])
  Y, X, _  = select_variables(df, [:T], x_names)
  fit = olsRegression(X ./10, Y, nothing, cluster)
  inference(fit)

  #Specification 3:
  x_names = vcat([:mean_grad_new],controlVariables[:controlVariableList])
  Y, X, _  = select_variables(df, [:T], x_names)
  fit = olsRegression(X  ./10 , Y, FixedEffects, cluster)
  fit.β
  inference(fit)
  
  #Specification 4:
  x_names = vcat([:mean_grad_new],controlVariables[:controlVariableList],controlVariables[:additionalControls])
  Y, X, _  = select_variables(df, [:T], x_names)
  fit = olsRegression(X ./10 , Y, FixedEffects, cluster)
  inference(fit)
  
  return table3
end






#----------------------------------------------------------------------------------------
# Tweak dataset to get exact point estimates:
#----------------------------------------------------------------------------------------


#Define some more rows first (column names / nObs row):
nCommunities = size(census_data,1)
nCommunitiesTreated = nrow(census_data[census_data[:,:T].==1,:])
nCommunitiesControl = nCommunities - nCommunitiesTreated

#Data modify:
#-----------

census_data = @chain census_data begin
    subset(:largeareas => ByRow(==(1)), skipmissing=true)
end
census_data[:,:kms_to_subs0] = census_data[:,:kms_to_subs0] ./10 ;
census_data[:,:baseline_hhdens0] = census_data[:,:baseline_hhdens0]./10;
census_data[:,:prop_indianwhite0] = census_data[:,:prop_indianwhite0]./10;
census_data[:,:kms_to_road0] = census_data[:,:kms_to_road0]./10 ;
census_data[:,:kms_to_town0] = census_data[:,:kms_to_town0]./10;


