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
using QuadGK
using ForwardDiff

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


# Part E: Estimate ATE, ATU, and ATT using MTE.
# - Restrict attention to parametric specifications of MTR
# - 1st estimating the treatment selection equation in (2) as a probit model to obtain estimates of the propensity score
# - 2nd modeling 𝐾(𝑝) as a polynomial in 𝑝 of degree k and estimating the outcome equation

# This one uses a probit built in function, planning to change it later: 

function get_propensity_score(df,covs)

    if length(covs) > 1
        reg = @formula(connected ~ treat_full + treat_low + treat_med + female + age_resp + educ + bank + housing + asset_value + energyspending + busia + market + transearly + connected_rate + population)
    else
        reg = @formula(connected ~ treat_full + treat_low + treat_med)    
    end

    probit = glm(reg, df, Binomial(), ProbitLink())

    propensity = predict(probit,df)

    return propensity

end


function get_mtr_parametrization(df, depvar=:hours_worked, covariates=:constant, treatment=:connected, K=1)

    propensity = get_propensity_score(df, [covariates])

    # Arrange all variables in matrix form:
    Y, X, D, _  = select_variables(df, [depvar], [covariates], [treatment])

    # Polinomial Approximation:
    propensityPoly = zeros(size(Y,1),K)
    for ii = 1:K propensityPoly[:,ii] = propensity.^ii end

    W = hcat(X, propensityPoly)
            
    # Parameters for MTR on identified segment:

    ols_d1 = olsRegression(Y[D[:].==1],W[D[:].==1,:],nothing,nothing,false)

    ols_d0 = olsRegression(Y[D[:].==0],W[D[:].==0,:],nothing,nothing,false)    

    β1 = ols_d1.β

    β0 = ols_d0.β
    
    return β1, β0, D, mean(propensity)
end



function discrete_mte_approximation_p1_nocovariates(df, depvar=:hours_worked, covariates=:constant, treatment=:connected, E_prte=0.424669)
    
    # Parametrize function as in Brinch Wiswall Mogstad (2017):
    β1, β0, D, E_p = get_mtr_parametrization(df, depvar, covariates, treatment, 1)

    k1(p) = β1[1] + β1[2]*p
    
    k0(p) = β0[1] + β0[2]*p

    # Now that we have all coefficients define function E[Y | P = p]:    
    mtr_approximation(p) = (k1(p))*p + (k0(p))*(1-p)

    mte_approximation(p) = ForwardDiff.derivative(mtr_approximation, p)


    # Weights:
    ω_ate(p)  = 1

    ω_att(p)  = ((1 - p) / mean(D))

    ω_atu(p)  = (p/ (1-mean(D)))

    # For ATE:
    mte_approximation_ate(p) = mte_approximation(p) * ω_ate(p)

    # For ATT:
    mte_approximation_att(p) = mte_approximation(p) * ω_att(p)

    # For ATU:
    mte_approximation_atu(p) = mte_approximation(p) * ω_atu(p)


    # Compute common estimates:
    ATE = mtr_approximation(1) - mtr_approximation(0) # This uses mtr

    ATE = quadgk(mte_approximation_ate, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    ATT = quadgk(mte_approximation_att, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    ATU = quadgk(mte_approximation_atu, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    PRTE = (mtr_approximation(E_prte) - mtr_approximation(E_p))/(E_prte-E_p) # This uses mtr

    # println("ATE: ",ATE, "\nATT: ",ATT,"\nATU: ",ATU)

    p_base = 0:0.01:1 ; MTR = mtr_approximation.(p_base)

    return ATE, ATT, ATU, PRTE, MTR, p_base 

end


function discrete_mte_approximation_p2_nocovariates(df, depvar=:hours_worked, covariates=:constant, treatment=:connected, E_prte=0.424669)
    
    # Parametrize function as in Brinch Wiswall Mogstad (2017):
    β1, β0, D, E_p = get_mtr_parametrization(df, depvar, covariates, treatment, 2)

    k1(p) = β1[1] + β1[2]*p + β1[3]*p^2
    
    k0(p) = β0[1] + β0[2]*p + β0[3]*p^2

    # Now that we have all coefficients define function E[Y | P = p]:    
    mtr_approximation(p) = (k1(p))*p + (k0(p))*(1-p)

    mte_approximation(p) = ForwardDiff.derivative(mtr_approximation, p)

    # Weights:
    ω_ate(p)  = 1

    ω_att(p)  = ((1 - p) / mean(D))

    ω_atu(p)  = (p/ (1-mean(D)))

    # For ATE:
    mte_approximation_ate(p) = mte_approximation(p) * ω_ate(p)

    # For ATT:
    mte_approximation_att(p) = mte_approximation(p) * ω_att(p)

    # For ATU:
    mte_approximation_atu(p) = mte_approximation(p) * ω_atu(p)


    # Compute common estimates:
    ATE = mtr_approximation(1) - mtr_approximation(0) # This uses mtr

    ATE = quadgk(mte_approximation_ate, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    ATT = quadgk(mte_approximation_att, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    ATU = quadgk(mte_approximation_atu, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    PRTE = (mtr_approximation(E_prte) - mtr_approximation(E_p))/(E_prte-E_p) # This uses mtr

    # println("ATE: ",ATE, "\nATT: ",ATT,"\nATU: ",ATU)

    p_base = 0:0.01:1 ; MTR = mtr_approximation.(p_base)

    return ATE, ATT, ATU, PRTE, MTR, p_base 

end


function discrete_mte_approximation_p3_nocovariates(df, depvar=:hours_worked, covariates=:constant, treatment=:connected, E_prte=0.424669)
    
    # Parametrize function as in Brinch Wiswall Mogstad (2017):
    β1, β0, D, E_p = get_mtr_parametrization(df, depvar, covariates, treatment, 3)

    k1(p) = β1[1] + β1[2]*p + β1[3]*p^2 + β1[4]*p^3
    
    k0(p) = β0[1] + β0[2]*p + β0[3]*p^2 + β0[4]*p^3

    # Now that we have all coefficients define function E[Y | P = p]:    
    mtr_approximation(p) = (k1(p))*p + (k0(p))*(1-p)

    mte_approximation(p) = ForwardDiff.derivative(mtr_approximation, p)

    # Weights:
    ω_ate(p)  = 1

    ω_att(p)  = ((1 - p) / mean(D))

    ω_atu(p)  = (p/ (1-mean(D)))

    # For ATE:
    mte_approximation_ate(p) = mte_approximation(p) * ω_ate(p)

    # For ATT:
    mte_approximation_att(p) = mte_approximation(p) * ω_att(p)

    # For ATU:
    mte_approximation_atu(p) = mte_approximation(p) * ω_atu(p)

    # Compute common estimates:
    ATE = mtr_approximation(1) - mtr_approximation(0) # This uses mtr

    ATE = quadgk(mte_approximation_ate, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    ATT = quadgk(mte_approximation_att, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    ATU = quadgk(mte_approximation_atu, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    PRTE = (mtr_approximation(E_prte) - mtr_approximation(E_p))/(E_prte-E_p) # This uses mtr

    # println("ATE: ",ATE, "\nATT: ",ATT,"\nATU: ",ATU)

    p_base = 0:0.01:1 ; MTR = mtr_approximation.(p_base)

    return ATE, ATT, ATU, PRTE, MTR, p_base 

end


function discrete_mte_approximation_p4_nocovariates(df, depvar=:hours_worked, covariates=:constant, treatment=:connected, E_prte=0.424669)
    
    # Parametrize function as in Brinch Wiswall Mogstad (2017):
    β1, β0, D, E_p = get_mtr_parametrization(df, depvar, covariates, treatment, 4)

    k1(p) = β1[1] + β1[2]*p + β1[3]*p^2 + β1[4]*p^3 + β1[5]*p^4
    
    k0(p) = β0[1] + β0[2]*p + β0[3]*p^2 + β0[4]*p^3 + β0[5]*p^4

    # Now that we have all coefficients define function E[Y | P = p]:    
    mtr_approximation(p) = (k1(p))*p + (k0(p))*(1-p)

    mte_approximation(p) = ForwardDiff.derivative(mtr_approximation, p)

    # Weights:
    ω_ate(p)  = 1

    ω_att(p)  = ((1 - p) / mean(D))

    ω_atu(p)  = (p/ (1-mean(D)))

    # ω_prte(p) = ((p-E_prte)/ (E_prte-E_p))

    # For ATE:
    mte_approximation_ate(p) = mte_approximation(p) * ω_ate(p)

    # For ATT:
    mte_approximation_att(p) = mte_approximation(p) * ω_att(p)

    # For ATU:
    mte_approximation_atu(p) = mte_approximation(p) * ω_atu(p)

    # For PRTE:
    # mte_approximation_prte(p) = mte_approximation(p) * ω_prte(p)


    # Compute common estimates:
    ATE = mtr_approximation(1) - mtr_approximation(0) # This uses mtr

    ATE = quadgk(mte_approximation_ate, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    ATT = quadgk(mte_approximation_att, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    ATU = quadgk(mte_approximation_atu, 0, 1, rtol=1e-8)[1] # This uses numerical integration

    # PRTE = quadgk(mte_approximation_prte, 0, 1, rtol=1e-8)[1]

    PRTE = (mtr_approximation(E_prte) - mtr_approximation(E_p))/(E_prte-E_p) # This uses mtr

    # println("ATE: ",ATE, "\nATT: ",ATT,"\nATU: ",ATU)

    p_base = 0:0.01:1 ; MTR = mtr_approximation.(p_base)

    return ATE, ATT, ATU, PRTE, MTR, p_base 

end


function bootstrap_se_model(df, depvar, discrete_mte_approximation)
    
    ATE_list = []; ATT_list = []; ATU_list = []; PRTE_list = []
    for _ in 1:1000
        index_b = rand(1:size(df,1),3000)
        df_b = df[index_b,:]
        try
            ATE, ATT, ATU, PRTE = discrete_mte_approximation(df_b, depvar)
            append!(ATE_list, ATE)
            append!(ATT_list, ATT)            
            append!(ATU_list, ATU)
            append!(PRTE_list, PRTE)
        catch
            println("Bootstrap error")
        end
    end

    ate_se = std(ATE_list)/sqrt(length(ATE_list)); 
    att_se = std(ATT_list)/sqrt(length(ATE_list))
    atu_se = std(ATU_list)/sqrt(length(ATE_list)); 
    prte_se = std(PRTE_list)/sqrt(length(PRTE_list)); 
    
    return ate_se, att_se, atu_se, prte_se

end




############################
# Implementation:
############################

# Part E: Estimate ATE, ATU, and ATT using MTE.
# - Restrict attention to parametric specifications of MTR
# - 1st estimating the treatment selection equation in (2) as a probit model to obtain estimates of the propensity score
# - 2nd modeling 𝐾(𝑝) as a polynomial in 𝑝 of degree k and estimating the outcome equation

df = load_dataset()
model_variables = define_variables()
covariate_names = [:constant]



# E) No covariates:
#------------------

# DEPVAR: Hours worked ... 

TreatmentEffects = zeros(4,8)  #[Order as follows: (ATE, ATT, ATU, PRTE)]

depvar = model_variables[:depvar_c3]

# Compile Treatment Effects:
TreatmentEffects[1,1], TreatmentEffects[2,1], TreatmentEffects[3,1], TreatmentEffects[4,1], _= discrete_mte_approximation_p1_nocovariates(df,depvar);
TreatmentEffects[1,3], TreatmentEffects[2,3], TreatmentEffects[3,3], TreatmentEffects[4,3], _= discrete_mte_approximation_p2_nocovariates(df,depvar);
TreatmentEffects[1,5], TreatmentEffects[2,5], TreatmentEffects[3,5], TreatmentEffects[4,5], _= discrete_mte_approximation_p3_nocovariates(df,depvar);
TreatmentEffects[1,7], TreatmentEffects[2,7], TreatmentEffects[3,7], TreatmentEffects[4,7], _= discrete_mte_approximation_p4_nocovariates(df,depvar);

# Compile Bootstrapped Standard Errors:
TreatmentEffects[1,2], TreatmentEffects[2,2], TreatmentEffects[3,2], TreatmentEffects[4,2]= bootstrap_se_model(df, depvar, discrete_mte_approximation_p1_nocovariates);
TreatmentEffects[1,4], TreatmentEffects[2,4], TreatmentEffects[3,4], TreatmentEffects[4,4]= bootstrap_se_model(df, depvar, discrete_mte_approximation_p2_nocovariates);
TreatmentEffects[1,6], TreatmentEffects[2,6], TreatmentEffects[3,6], TreatmentEffects[4,6]= bootstrap_se_model(df, depvar, discrete_mte_approximation_p3_nocovariates);
TreatmentEffects[1,8], TreatmentEffects[2,8], TreatmentEffects[3,8], TreatmentEffects[4,8]= bootstrap_se_model(df, depvar, discrete_mte_approximation_p4_nocovariates);


CSV.write("Q1_PE_hours_woked.csv",  Tables.table(round.(TreatmentEffects,sigdigits = 3)))



# DEPVAR: Life Satisfaction...

TreatmentEffects = zeros(4,8)  #[Order as follows: (ATE, ATT, ATU)]

depvar = model_variables[:depvar_d2]

# Compile Treatment Effects:
TreatmentEffects[1,1], TreatmentEffects[2,1], TreatmentEffects[3,1], TreatmentEffects[4,1], _= discrete_mte_approximation_p1_nocovariates(df,depvar);
TreatmentEffects[1,3], TreatmentEffects[2,3], TreatmentEffects[3,3], TreatmentEffects[4,3], _= discrete_mte_approximation_p2_nocovariates(df,depvar);
TreatmentEffects[1,5], TreatmentEffects[2,5], TreatmentEffects[3,5], TreatmentEffects[4,5], _= discrete_mte_approximation_p3_nocovariates(df,depvar);
TreatmentEffects[1,7], TreatmentEffects[2,7], TreatmentEffects[3,7], TreatmentEffects[4,7], _= discrete_mte_approximation_p4_nocovariates(df,depvar);

# Compile Bootstrapped Standard Errors:
TreatmentEffects[1,2], TreatmentEffects[2,2], TreatmentEffects[3,2], TreatmentEffects[4,2]= bootstrap_se_model(df, depvar, discrete_mte_approximation_p1_nocovariates);
TreatmentEffects[1,4], TreatmentEffects[2,4], TreatmentEffects[3,4], TreatmentEffects[4,4]= bootstrap_se_model(df, depvar, discrete_mte_approximation_p2_nocovariates);
TreatmentEffects[1,6], TreatmentEffects[2,6], TreatmentEffects[3,6], TreatmentEffects[4,6]= bootstrap_se_model(df, depvar, discrete_mte_approximation_p3_nocovariates);
TreatmentEffects[1,8], TreatmentEffects[2,8], TreatmentEffects[3,8], TreatmentEffects[4,8]= bootstrap_se_model(df, depvar, discrete_mte_approximation_p4_nocovariates);


CSV.write("Q1_PE_life_satisfaction.csv",  Tables.table(round.(TreatmentEffects,sigdigits = 3)))



# F) With covariates:
#--------------------
K = 2
covariates=[]
append!(covariates,model_variables[:household])
append!(covariates,model_variables[:comunity])
covariate_names = vcat(covariates, [:constant])
support  = (0.01, 1)

#Drop missings:
df_with_covs = df[:,vcat([model_variables[:depvar_c3], model_variables[:depvar_d2]], 
                            [model_variables[:treatment]], 
                            covariate_names,
                            model_variables[:instruments])]

df_with_covs = df_with_covs[all.(!ismissing, eachrow(df_with_covs)), :]


# G) Now consider policy where effective connection price is changed to $200
# - Construct estimate of per person policy relevant treatment relative to status quo
# - use two outcomes in (a) and parametric, point identified specifications of mtr as in e.

# Small tweak, instead of dichotomous variables I'll transform it in effective connection price:
# - MTE will remain the same, we will just change propensity score and compute PRTE

# Prte Experiment
probit_prte = glm(@formula(connected ~ effective_connection_price), df, Binomial(), ProbitLink())
df_pol = copy(df)
df_pol[:,:effective_connection_price] .= 200
E_prte = mean(predict(probit_prte,df_pol))

#Note: I defaulted E_prte value in the previous functions to automatically get the PRTE in all polinomial specifications
# That's why they appear in the tables








