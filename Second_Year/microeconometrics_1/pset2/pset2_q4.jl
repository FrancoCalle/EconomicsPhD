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
# Pkg.instantiate()
Pkg.add("Missings")
Pkg.add("FixedEffectModels")


#Load packages ...

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


function jacknife_iv_estimation(y, d, z, x=nothing, fixed_effects=nothing, intercept=false)

  # Make matrix for both first and second stage 
  n = size(y,1)
        
  if ~isnothing(x)
      z = hcat(z,x)
      d = hcat(d,x)
  end

  if ~isnothing(fixed_effects)
      z = hcat(z,fixed_effects)
      d = hcat(d,fixed_effects)
  end

  if intercept == true
      z = hcat(z,ones(n,1))
      d = hcat(d,ones(n,1))
  end

  k = size(d,2)

  # Operations
  ZᵀZ = z' * z
  L_i = diag(z * inv(ZᵀZ) * z')
  Π_hat = d' * z * inv(ZᵀZ)
  d_jive = (1 ./ (1 .- L_i)).*(z * Π_hat' .- L_i .* d)

  # Coefficient
  d_jiveᵀd = d_jive' * d
  d_jiveᵀY = d_jive' * y
  beta = inv(d_jiveᵀd) * d_jiveᵀY

  # Get residuals
  # res = Y - d * beta
  # U2 = res * res'
  # Sigma = inv(d_jiveᵀd) * (d_jive' * diag(U2) * d_jive) * inv(d_jiveᵀd)' * n/(n-k)
  # Sigma = sum(U2)/(n-k) * inv(d_jiveᵀd) * (d_jive' * d_jive) * inv(d_jiveᵀd)'

  return beta

end


function solve_part_c(depvar,df=census_data)

  instrument = [:mean_grad_new]
  treatment = [:T]
  FixedEffects = Matrix(reduce(hcat, [df[:,:dccode0].== fe for fe in unique(df[:,:dccode0])]))
  x_names = controlVariables[:controlVariableList]

  # Specification 5:
  Y, _, D, Z  = select_variables(census_data, depvar, x_names, treatment, instrument)
  β_5 = jacknife_iv_estimation(Y, D, Z, nothing, nothing, true)

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

  return Dict(:model5 => β_5, :model6 => β_6, :model7 => β_7, :model8 => β_8)

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

# Part A: Replicate tables
#--------------------------------

t3_results = generateTable3(census_data, controlVariables)
fem_t4_results = generateTable4(census_data, [:d_prop_emp_f], controlVariables)
male_t5_results = generateTable4(census_data, [:d_prop_emp_m], controlVariables)


npanels = length(male_t5_results)
r = length(male_t5_results[:panel_8].β)
table_results = zeros(r,npanels)
for ii in 1:8
  model = Symbol("panel_", ii)
  table_results[1:length(male_t5_results[model].β),ii] = male_t5_results[model].β
end

CSV.write("table4_male_results.csv",  Tables.table(round.(table_results,digits = 3)), writeheader=false)



# Part B: Anderson Rubin test finer grid
#---------------------------------------

# Dinkelman Grid: [-0.6: 0.05 : 1]
β_set_dinkelman = anderson_rubin_test(census_data, [:d_prop_emp_f], controlVariables, -.6, 1, 0.05)

# Finer Grid: [-0.6: 0.05 : 1]
β_set_finer = anderson_rubin_test(census_data, [:d_prop_emp_f], controlVariables, -.3, 0.6, 0.005)


# Part C: Use Jacknife
#---------------------------------------

results_female = solve_part_c([:d_prop_emp_f],census_data)
results_male = solve_part_c([:d_prop_emp_m], census_data)



r = size(results_female[:model8],1)
table_results = zeros(r,4)
for ii in 1:4
  model = Symbol("model", ii+4)
  table_results[1:length(results_female[model]),ii] = results_female[model]
end  
CSV.write("jive_female_results.csv",  Tables.table(round.(table_results,sigdigits = 3)), writeheader=false)


r = size(results_male[:model8],1)
table_results = zeros(r,4)
for ii in 1:4
  model = Symbol("model", ii+4)
  table_results[1:length(results_male[model]),ii] = results_male[model]
end  
CSV.write("jive_male_results.csv",  Tables.table(round.(table_results,sigdigits = 3)), writeheader=false)




