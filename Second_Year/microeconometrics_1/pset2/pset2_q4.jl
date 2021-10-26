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
include(joinpath(@__DIR__, "..", "..", "..","fc_toolkit.jl"))


function setVariableNames()
  
  controlVariableList = [:base_hhpovrate0,
                          :prop_head_f_a0, 
                          :sexratio0_a, 
                          :prop_indianwhite0, 
                          :kms_to_road0, 
                          :kms_to_town0,  
                          :prop_matric_m0,
                          :prop_matric_f0, 
                          :baseline_hhdens0, 
                          :kms_to_subs0]

  additionalControls = [:d_prop_waterclose, 
                        :d_prop_flush]

  controlVariables = Dict(:controlVariableList => controlVariableList, 
                          :additionalControls => additionalControls)

  return controlVariables

end


@doc """
Function that generates table 3 (First Stage):

df: is a dataframe which contains the variables used in regressions for table3
controlVariables: is a list which contains lists of names and tag of variables  
"""->
function generateTable3(df, controlVariables)

  cluster = CategoricalArray(df[:, :placecode0])
  FixedEffects = Matrix(reduce(hcat, [df[:,:dccode0].==fe for fe in unique(df[:,:dccode0])]))

  #Specification 1:
  Y, X, _  = select_variables(df, [:T], [:mean_grad_new])
  ols1 = olsRegression(Y, X ./10 ,  nothing, cluster)
  ols1 = inference(ols1)

  #Specification 2:
  x_names = vcat([:mean_grad_new],controlVariables[:controlVariableList])
  Y, X, _  = select_variables(df, [:T], x_names)
  ols2 = olsRegression(Y, X ./10,  nothing, cluster)
  ols2 = inference(ols2)

  #Specification 3:
  x_names = vcat([:mean_grad_new],controlVariables[:controlVariableList])
  Y, X, _  = select_variables(df, [:T], x_names)
  ols3 = olsRegression(Y, X  ./10 , FixedEffects, cluster)
  ols3 = inference(ols3)
  
  #Specification 4:
  x_names = vcat([:mean_grad_new],controlVariables[:controlVariableList],controlVariables[:additionalControls])
  Y, X, _  = select_variables(df, [:T], x_names)
  ols4 = olsRegression(Y, X ./10 , FixedEffects, cluster)
  ols4 = inference(ols4)
  
  #Compile all in dictionary:
  results = Dict( :panel_1 => ols1,
                  :panel_2 => ols2,
                  :panel_3 => ols3,
                  :panel_4 => ols4)

  return results

end


@doc """
  Function that generates table 3 (First Stage):

  df: is a dataframe which contains the variables used in regressions for table3
  depvar: a 1x1 list which contains the dependent variable
  controlVariables: is a list which contains lists of names  
"""->
function generateTable4(df, depvar, controlVariables)
  # TODO: work on inference function for tsls ...

  instrument = [:mean_grad_new]
  treatment = [:T]
  cluster = CategoricalArray(df[:, :placecode0])
  FixedEffects = Matrix(reduce(hcat, [df[:,:dccode0].==fe for fe in unique(df[:,:dccode0])]))


  #Ordinary Least Squares [Panels 1 - 4] CHEFF KISS: (for point estimates)
  #-----------------------------------------------------------------------
  #Specification 1:
  x_names = treatment
  Y, X, _, _  = select_variables(df, depvar, x_names)
  ols1 = olsRegression(Y, X, nothing, cluster)
  # inference(fit)

  #Specification 2:
  x_names = vcat(treatment, controlVariables[:controlVariableList])
  Y, X, _, _  = select_variables(df, depvar, x_names)
  ols2 = olsRegression(Y, X, nothing, cluster)
  # inference(fit)

  #Specification 3:
  x_names = vcat(treatment, controlVariables[:controlVariableList])
  Y, X, _, _  = select_variables(df, depvar, x_names)
  ols3 = olsRegression(Y, X, FixedEffects, cluster)
  # inference(fit)
  
  #Specification 4:
  x_names = vcat(treatment, controlVariables[:controlVariableList], controlVariables[:additionalControls])
  Y, X, _, _  = select_variables(df, depvar, x_names)
  ols4 = olsRegression(Y, X, FixedEffects, cluster)

  #Two stage least squares [Panels 5 - 8] CHEFF KISS: (for point estimates)
  #------------------------------------------------------------------------

  #Specification 1:
  Y, _, D, Z  = select_variables(df, depvar, nothing, treatment, instrument)
  tsls1 = tsls_regression(Y, D, Z)
  # inference(fit)

  #Specification 2:
  x_names = controlVariables[:controlVariableList]
  Y, X, D, Z  = select_variables(df, depvar, x_names, treatment, instrument)
  tsls2 = tsls_regression(Y, D, Z, X)
  # inference(fit)

  #Specification 3:
  x_names = controlVariables[:controlVariableList]
  Y, X, D, Z  = select_variables(df, depvar, x_names, treatment, instrument)
  tsls3 = tsls_regression(Y, D, Z, X, FixedEffects)
  # inference(fit)
  
  #Specification 4:
  x_names = vcat(controlVariables[:controlVariableList],controlVariables[:additionalControls])
  Y, X, D, Z  = select_variables(df, depvar, x_names, treatment, instrument)
  tsls4 = tsls_regression(Y, D, Z, X, FixedEffects)
  # inference(fit)
  
  results = Dict( :panel_1 => ols1,
                  :panel_2 => ols2,
                  :panel_3 => ols3,
                  :panel_4 => ols4,
                  :panel_5 => tsls1,
                  :panel_6 => tsls2,
                  :panel_7 => tsls3,
                  :panel_8 => tsls4 )

  return results

end


function anderson_rubin_test(df, depvar, controlVariables, min_β=-.3, max_β=0.6, δ =0.05)

  # Define symbols we will use...
  instrument = [:mean_grad_new]
  treatment = [:T]
  cluster = CategoricalArray(df[:, :placecode0])
  FixedEffects = Matrix(reduce(hcat, [df[:,:dccode0].==fe for fe in unique(df[:,:dccode0])]))
  
  # Get data in arrays...
  x_names = vcat(treatment, controlVariables[:controlVariableList], controlVariables[:additionalControls])
  Y, X, _, Z  = select_variables(df, depvar, x_names, nothing, instrument)

  # Define grid search...
  β_list = Array(min_β:δ:max_β)
  
  # Compute Anderson Rubin test for each element in grid...
  function ar_test(b, y=Y, x=X, z = Z, FixedEffects=FixedEffects, cluster=cluster)
    resid = y - b .* x[:,1]
    artest_m = olsRegression(resid, hcat(z,x[:,2:end]), FixedEffects, cluster)
    artest_m = inference(artest_m)
    return artest_m.p[1]
  end

  # Obtain p values...
  p_list = ar_test.(β_list)

  # Index betas that reject null at 0.05
  return β_list[p_list.>=0.05]

end



function jacknife_iv_estimation(Y, D, Z, X=nothing, fe=nothing, constant=nothing)

  # Make matrix for both first and second stage 
  n = size(Y,1)

  if ~isnothing(X)
    Z = hcat(Z,X)
    X = hcat(D,X)
  end

  if ~isnothing(fe)
    constant = false
    Z = hcat(Z,fe)
    X = hcat(X,fe)
  end

  if constant == true
    X = hcat(X,ones(n,1))
    Z = hcat(Z,ones(n,1))
  end

  k = size(X,2)

  # Operations
  ZᵀZ = Z' * Z
  L_i = diag(Z * inv(ZᵀZ) * Z')
  Π_hat = X' * Z * inv(ZᵀZ)
  X_jive = (1 ./ (1 .- L_i)).*(Z * Π_hat' .- L_i .* X)

  # Coefficient
  X_jiveᵀX = X_jive' * X
  X_jiveᵀY = X_jive' * Y
  beta = inv(X_jiveᵀX) * X_jiveᵀY

  # Get residuals
  # res = Y - X * beta
  # U2 = res * res'

  # if varcov == "het"
  #   Sigma = inv(X_jiveᵀX) * (X_jive' * diag(U2) * X_jive) * inv(X_jiveᵀX)' * n/(n-k)
  # elseif varcov == "hom"  
  #   Sigma = sum(U2)/(n-k) * inv(X_jiveᵀX) * (X_jive' * X_jive) * inv(X_jiveᵀX)'
  # end

  return beta

end

# Small data clean and execute the code:
#---------------------------------------

census_data = DataFrame(load("matched_censusdata.dta"));

census_data = @chain census_data begin
  subset(:largeareas => ByRow(==(1)), skipmissing=true)
end
census_data[:,:kms_to_subs0] = census_data[:,:kms_to_subs0] ./10 ;
census_data[:,:baseline_hhdens0] = census_data[:,:baseline_hhdens0]./10;
census_data[:,:prop_indianwhite0] = census_data[:,:prop_indianwhite0]./10;
census_data[:,:kms_to_road0] = census_data[:,:kms_to_road0]./10 ;
census_data[:,:kms_to_town0] = census_data[:,:kms_to_town0]./10;

controlVariables  = setVariableNames()
t3_results = generateTable3(census_data, controlVariables)
fem_t4_results = generateTable4(census_data, [:d_prop_emp_f], controlVariables)
fem_t5_results = generateTable4(census_data, [:d_prop_emp_m], controlVariables)


# Dinkelman Grid: [-0.6: 0.05 : 1]
β_set_dinkelman = anderson_rubin_test(census_data, [:d_prop_emp_f], controlVariables, -.6, 1, 0.05)

# Finer Grid: [-0.6: 0.05 : 1]
β_set_finer = anderson_rubin_test(census_data, [:d_prop_emp_f], controlVariables, -.3, 0.6, 0.005)


# Part C: Use Jacknife
#---------------------------------------

function solve_part_c(depvar)

  instrument = [:mean_grad_new]
  treatment = [:T]
  FixedEffects = Matrix(reduce(hcat, [df[:,:dccode0].==fe for fe in unique(df[:,:dccode0])]))
  x_names = controlVariables[:controlVariableList]

  # Specification 6:
  Y, X, D, Z  = select_variables(census_data, depvar, x_names, treatment, instrument)
  β_6 = jacknife_iv_estimation(Y, D, Z, X, nothing, true)

  # Specification 7:
  Y, X, D, Z  = select_variables(census_data, depvar, x_names, treatment, instrument)
  β_7 = jacknife_iv_estimation(Y, D, Z, X, FixedEffects)

  # Specification 8:
  x_names = vcat(controlVariables[:controlVariableList],controlVariables[:additionalControls])
  Y, X, D, Z  = select_variables(census_data, depvar, x_names, treatment, instrument)
  β_8 = jacknife_iv_estimation(Y, D, Z, X, FixedEffects)

  return Dict(:model6 => β_6, :model7 => β_7, :model8 => β_8)

end


results_female = solve_part_c([:d_prop_emp_f])
results_male = solve_part_c([:d_prop_emp_m])