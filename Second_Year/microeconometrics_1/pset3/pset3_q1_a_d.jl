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
using GLM

using JuMP, GLPK # for linear optimization

include(joinpath(@__DIR__,"..", "..", "..","fc_toolkit.jl"))

function load_dataset()

    df = DataFrame(CSV.File("lmw2019jpe-clean.csv"))

    df[:, :effective_connection_price] .= 398
    df[df.treat_full.==1, :effective_connection_price] .= 0
    df[df.treat_med.==1, :effective_connection_price] .= 171
    df[df.treat_low.==1, :effective_connection_price] .= 284
    
    df[:,:constant] .= 1

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

    household = [:female, :age_resp, :educ, :bank, :housing, :asset_value, :energyspending, :round2]

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


function compute_results_part_a()
    
    # Replicate C3: Total Hours Worked :
    #-----------------------------------
    println("Total hours worked")
    depvar = [model_variables[:depvar_c3]]

    println("ITT:")
    Y, X, _, CL= select_variables(df, depvar, x_vars, nothing, model_variables[:cluster])
    CL = CategoricalArray(CL)
    ols_c3 = olsRegression(Y,X, nothing, CL[:])
    ols_c3 = inference(ols_c3, "clust")

    println("TOT:")
    Y, X, D, Z  = select_variables(df, depvar, covariates, treatment, instrument)
    tsls_c3 = tsls_regression(Y, D, Z, X)
    tsls_c3 = inference(tsls_c3)

    # Replicate D2 TOT: Life Satisfaction:
    #-------------------------------------
    println("Life Satisfaction")
    depvar = [model_variables[:depvar_d2]]

    println("ITT:")
    Y, X, _, CL= select_variables(df, depvar, x_vars, nothing, model_variables[:cluster])
    CL = CategoricalArray(CL)
    ols_d2 = olsRegression(Y,X, nothing, CL[:])
    ols_d2 = inference(ols_d2, "clust")

    println("TOT:")
    Y, X, D, Z  = select_variables(df, depvar, covariates, treatment, instrument)
    tsls_d2 = tsls_regression(Y, D, Z, X)
    tsls_d2 = inference(tsls_d2)

    # Compile results:
    results = zeros(4,2)
    #ITT
    results[1,1] = ols_c3.β[1]
    results[2,1] = ols_c3.se[1]

    results[3,1] = ols_d2.β[1]
    results[4,1] = ols_d2.se[1]

    #TOT:
    results[1,2] = tsls_c3.β[1]
    results[2,2] = tsls_c3.se[1]

    results[3,2] = tsls_d2.β[1]
    results[4,2] = tsls_d2.se[1]

    return results, ols_c3, tsls_c3, ols_d2, tsls_d2
    
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
    ols_c3 = olsRegression(Y,X, nothing, CL[:])
    ols_c3 = inference(ols_c3, "clust")

    println("TOT:")
    Y, _, D, Z  = select_variables(df, depvar, [:constant], treatment, instrument)
    tsls_c3 = tsls_regression(Y, D, Z)
    tsls_c3 = inference(tsls_c3)

    println(mean(Y[X[:].==0]), " ",std(Y[X[:].==0]))

    # Replicate D2 TOT: Life Satisfaction:
    #-------------------------------------
    println("Life Satisfaction")
    depvar = [model_variables[:depvar_d2]]

    println("ITT:")
    Y, X, _, CL= select_variables(df, depvar, [itt_treatment], nothing, model_variables[:cluster])
    CL = CategoricalArray(CL)
    ols_d2 = olsRegression(Y,X, nothing, CL[:])
    ols_d2 = inference(ols_d2, "clust")

    println("TOT:")
    Y, _, D, Z  = select_variables(df, depvar, [:constant], treatment, instrument)
    tsls_d2 = tsls_regression(Y, D, Z)
    tsls_d2 = inference(tsls_d2)

    println(mean(Y[X.==0]),std(Y[X.==0]))

    # Compile results:
    results = zeros(4,2)
    #ITT
    results[1,1] = ols_c3.β[1]
    results[2,1] = ols_c3.se[1]

    results[3,1] = ols_d2.β[1]
    results[4,1] = ols_d2.se[1]

    #TOT:
    results[1,2] = tsls_c3.β[1]
    results[2,2] = tsls_c3.se[1]

    results[3,2] = tsls_d2.β[1]
    results[4,2] = tsls_d2.se[1]

    return results, ols_c3, tsls_c3, ols_d2, tsls_d2

end


# Parts A and D:
#---------------

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

results_a,_ = compute_results_part_a()
results_d,_ = compute_results_part_d()

CSV.write("Q1_PA_ITT_TOT.csv",  Tables.table(round.(results_a,sigdigits = 3)))
CSV.write("Q1_PD_ITT_TOT.csv",  Tables.table(round.(results_d,sigdigits = 3)))













