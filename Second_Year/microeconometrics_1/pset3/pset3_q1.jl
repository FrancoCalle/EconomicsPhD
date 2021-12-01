using Pkg
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

function load_dataset()

    df = DataFrame(CSV.File("lmw2019jpe-clean.csv"))

    return df 

end


function define_variables()
    depvar_c3 = :hours_worked 
    depvar_d2 = :life_index
    instruments = [:treat_low, :treat_med, :treat_full] 
    treatment = :connected
    itt_treatment = :treat_full

    # Covariates:
    covariates_comunity = [:busia, :market, :transearly, :connected_rate, :population]
    household = [:female, :base_age, :educ, :bank, :housing, :asset_value, :energyspending]

    cluster = [:siteno]

    return Dict(:depvar_c3 => depvar_c3, 
                :depvar_d2 => depvar_d2, 
                :instruments => instruments,
                :treatment => treatment,
                :itt_treatment => itt_treatment,
                :comunity => covariates_comunity,
                :household => household,
                :cluster => cluster)
end

df = load_dataset()

model_variables = define_variables()

# Part A: Replicate cols (1) and (3) rows C3 and D2 of table 3 

# Specification for ITT only uses "treat_full", meanwhile TOT instruments electrification status vs all subsidies
instrument = model_variables[:instruments]
treatment = [model_variables[:treatment]]
covariates=[]
append!(covariates,model_variables[:household])
append!(covariates,model_variables[:comunity])
itt_treatment = model_variables[:itt_treatment]

x_vars = []
x_vars = append!(x_vars, [itt_treatment])
x_vars = append!(x_vars, covariates)


function compute_results_part_a()
    
    # Replicate C3: Total Hours Worked :
    #-----------------------------------
    println("Total hours worked")
    depvar = [model_variables[:depvar_c3]]

    println("ITT:")
    Y, X, _, CL= select_variables(df, depvar, x_vars, nothing, model_variables[:cluster])
    CL = CategoricalArray(CL)
    ols1 = olsRegression(Y,X, nothing, CL[:])
    ols1 = inference(ols1)
    ols1.se
    ols1.β


    println("TOT:")
    Y, X, D, Z  = select_variables(df, depvar, covariates, treatment, instrument)
    tsls1 = tsls_regression(Y, D, Z, X)
    tsls1.β

    # Replicate D2 TOT: Life Satisfaction:
    #-------------------------------------
    println("Life Satisfaction")
    depvar = [model_variables[:depvar_d2]]

    println("ITT:")
    Y, X, _, CL= select_variables(df, depvar, x_vars, nothing, model_variables[:cluster])
    CL = CategoricalArray(CL)
    ols1 = olsRegression(Y,X, nothing, CL[:])
    ols1.β
    inference(ols1)

    println("TOT:")
    Y, X, D, Z  = select_variables(df, depvar, covariates, treatment, instrument)
    tsls1 = tsls_regression(Y, D, Z, X)
    tsls1.β
    
end



# Part D: Show that results in (a) do not change much by not controlling for covariates
function compute_results_part_d()
    
    # Replicate C3 TOT: Life Satisfaction:
    #-------------------------------------
    println("Total hours worked")
    depvar = [model_variables[:depvar_c3]]

    println("ITT:")
    Y, X, _, CL= select_variables(df, depvar, [itt_treatment], nothing, model_variables[:cluster])
    CL = CategoricalArray(CL)
    ols1 = olsRegression(Y,X, nothing, CL[:])
    ols1 = inference(ols1)
    ols1.β
    ols1.se


    println("TOT:")
    Y, X, D, Z  = select_variables(df, depvar, covariates, treatment, instrument)
    tsls1 = tsls_regression(Y, D, Z)
    tsls1.β

    # Replicate D2 TOT: Life Satisfaction:
    #-------------------------------------
    println("Life Satisfaction")
    depvar = [model_variables[:depvar_d2]]

    println("ITT:")
    Y, X, _, CL= select_variables(df, depvar, [itt_treatment], nothing, model_variables[:cluster])
    CL = CategoricalArray(CL)
    ols1 = olsRegression(Y,X, nothing, CL[:])
    ols1.β
    inference(ols1)

    println("TOT:")
    Y, X, D, Z  = select_variables(df, depvar, covariates, treatment, instrument)
    tsls1 = tsls_regression(Y, D, Z)
    tsls1.β
 
end





# Part E: Estimate ATE, ATU, and ATT using MTE.
# - Restrict attention to parametric specifications of MTR
# - Use bootstrap to produce SE.
# - Examine sensitivity of results to parametrization of MTR

# Part F: Repeat E while controlling for the same set of cocvs as in A.

# Part G: Consider policy where effective connection price is changed to 200 USD for all.
# - Construct estimate of per-person policy-relevant treatment effect of new policy
# - Use bootstrap to produce SE or CI.

# Part H: Pick prior lower and upper bounds for the two outcomes in A and  discuss choices.
# - Compute nonparametric bounds on the ATU that use all info in E[Y_i | D_i, Z_i].

# Parametrize MTR function using Bernstein polynomials, omiting covariates and recompute bounds. 

















