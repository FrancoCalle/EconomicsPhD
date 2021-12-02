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

include(joinpath(@__DIR__,"..", "..", "..","fc_toolkit.jl"))

function load_dataset()

    df = DataFrame(CSV.File("lmw2019jpe-clean.csv"))
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
    ols_c3 = olsRegression(Y,X, nothing, CL[:])
    ols_c3 = inference(ols1)
    ols_c3.se
    ols_c3.β


    println("TOT:")
    Y, X, D, Z  = select_variables(df, depvar, covariates, treatment, instrument)
    tsls_c3 = tsls_regression(Y, D, Z, X)
    tsls_c3.β

    # Replicate D2 TOT: Life Satisfaction:
    #-------------------------------------
    println("Life Satisfaction")
    depvar = [model_variables[:depvar_d2]]

    println("ITT:")
    Y, X, _, CL= select_variables(df, depvar, x_vars, nothing, model_variables[:cluster])
    CL = CategoricalArray(CL)
    ols_d2 = olsRegression(Y,X, nothing, CL[:])
    ols_d2.β
    inference(ols1)

    println("TOT:")
    Y, X, D, Z  = select_variables(df, depvar, covariates, treatment, instrument)
    tsls_d2 = tsls_regression(Y, D, Z, X)
    tsls_d2.β

    return ols_c3, tsls_c3, ols_d2, tsls_d2
    
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
    ols_c3 = inference(ols1)
    ols_c3.β
    ols_c3.se

    println("TOT:")
    Y, X, D, Z  = select_variables(df, depvar, covariates, treatment, instrument)
    tsls_c3 = tsls_regression(Y, D, Z)
    tsls_c3.β

    # Replicate D2 TOT: Life Satisfaction:
    #-------------------------------------
    println("Life Satisfaction")
    depvar = [model_variables[:depvar_d2]]

    println("ITT:")
    Y, X, _, CL= select_variables(df, depvar, [itt_treatment], nothing, model_variables[:cluster])
    CL = CategoricalArray(CL)
    ols_d2 = olsRegression(Y,X, nothing, CL[:])
    ols_d2.β
    inference(ols1)

    println("TOT:")
    Y, X, D, Z  = select_variables(df, depvar, covariates, treatment, instrument)
    tsls_d2 = tsls_regression(Y, D, Z)
    tsls_d2.β
    
    return ols_c3, tsls_c3, ols_d2, tsls_d2

end





# Part E: Estimate ATE, ATU, and ATT using MTE.
# - Restrict attention to parametric specifications of MTR
# - 1st estimating the treatment selection equation in (2) as a probit model to obtain estimates of the propensity score
# - 2nd modeling 𝐾(𝑝) as a polynomial in 𝑝 of degree k and estimating the outcome equation

# This one uses a probit built in function, planning to change it later: 

function get_propensity_score()

    reg = @formula(connected ~ treat_full + treat_low + treat_med + female + base_age + educ + bank + housing + asset_value + energyspending + busia + market + transearly + connected_rate + population)

    probit = glm(reg, df, Binomial(), ProbitLink())

    propensity = predict(probit)

    return propensity

end


function marginalEffect(β10, α, XPr, polinomialApproximation, propensity)
    
    C = Array(2:K)

    dkdp = C'.*(polinomialApproximation./propensity)
    
    mte = XPr*β10 .+ dkdp*α

    return mte

end


function polinomial_approximation(K, propensity)

    polinomialApproximation = zeros(length(propensity), K-1)

    for kk in 1:K-1 polinomialApproximation[:,kk] = propensity.^(kk+1) end
    
    return polinomialApproximation

end


function unpack_parameters(parameters, nCovs)

    β0 = parameters[1:nCovs]

    β10 = parameters[nCovs+1:2*nCovs]

    α  = parameters[2*nCovs+1 : end]
    
    return β0, β10, α

end


function get_marginal_response(depvar, covariates, K)

    propensity = get_propensity_score()

    polinomialApproximation = polinomial_approximation(K, propensity)
    
    # Arrange all variables in matrix form:
    Y, X, D, _  = select_variables(df, [depvar], covariates, [:connected])

    # Recall function to estimate: Y = X β_0 + X(β_1 - β_0)̂p + K(̂p) + ϵ
    
    XX = X  # # nX = size(X,1); XX = hcat(X, ones(nX,1))

    XPr = XX .* propensity

    W = hcat(XX, XPr, polinomialApproximation)

    ols1 = olsRegression(Y,W,nothing,nothing,false)

    _, β10, α = unpack_parameters(ols1.β, size(XX,2))

    mte = marginalEffect(β10, α, XPr, polinomialApproximation, propensity)

    return mte, propensity, Y, XX, D

end


function get_average_treatmet(mte, D)

    ATE = mean(mte)
    ATU = mean(mte[D[:] .== 0])
    ATT = mean(mte[D[:] .== 1])

    return ATE, ATU, ATT

end

K = 3
depvar = model_variables[:depvar_c3]
covariate_names = vcat(covariates, [:constant])
mte, propensity, XX, Y, D = get_marginal_response(depvar, covariate_names, K)


mean(mte)
histogram(propensity)
histogram(mte)

mean(mte)

# And derivative is: ∂E[Y | X=x, P(Z) = p]/∂p = X(β1 - β0) + ∂ K(p) / ∂ p

# histogram(predict(probit))




# - Use bootstrap to produce SE.

# - Examine sensitivity of results to parametrization of MTR

# Part F: Repeat E while controlling for the same set of cocvs as in A.

# Part G: Consider policy where effective connection price is changed to 200 USD for all.
# - Construct estimate of per-person policy-relevant treatment effect of new policy
# - Use bootstrap to produce SE or CI.

# Part H: Pick prior lower and upper bounds for the two outcomes in A and  discuss choices.
# - Compute nonparametric bounds on the ATU that use all info in E[Y_i | D_i, Z_i].

# Parametrize MTR function using Bernstein polynomials, omiting covariates and recompute bounds. 

















