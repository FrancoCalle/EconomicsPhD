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


function marginalEffect(β10, α, XX, polinomialApproximation, propensity)

    #  Derivative is: ∂E[Y | X=x, P(Z) = p]/∂p = X(β1 - β0) + ∂ K(p) / ∂ p

    C = Array(2:K)

    dkdp = C'.*(polinomialApproximation./propensity)
    
    mte = XX*β10 .+ dkdp*α

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


function get_marginal_response(df, depvar, covariates, K, h_lo, h_up)

    propensity = get_propensity_score(df, covariates)

    polinomialApproximation = polinomial_approximation(K, propensity)
    
    # Arrange all variables in matrix form:
    Y, X, D, _  = select_variables(df, [depvar], covariates, [:connected])

    get_indentified_segment(X, propensity, h_lo, h_up) = X[(propensity .> h_lo) .& (propensity .< h_up),:]

    # Recall function to estimate: Y = X β_0 + X(β_1 - β_0)̂p + K(̂p) + ϵ
    XX = X  # # nX = size(X,1); XX = hcat(X, ones(nX,1))

    XPr = XX .* propensity

    W = hcat(XX, XPr, polinomialApproximation)

    # Get identified segment for all variables:
    
    Wπ = get_indentified_segment(W, propensity, h_lo, h_up)
    
    XXπ = get_indentified_segment(XX, propensity, h_lo, h_up)

    Yπ = get_indentified_segment(Y, propensity, h_lo, h_up)

    Dπ = get_indentified_segment(D, propensity, h_lo, h_up)

    XPrπ = get_indentified_segment(XPr, propensity, h_lo, h_up)
    
    propensityπ = get_indentified_segment(propensity, propensity, h_lo, h_up)

    polinomialApproximationπ = get_indentified_segment(polinomialApproximation, propensity, h_lo, h_up)
        
    # Compute mte on identified segment:

    ols1 = olsRegression(Yπ,Wπ,nothing,nothing,false)

    _, β10, α = unpack_parameters(ols1.β, size(XX,2))

    mte = marginalEffect(β10, α, XXπ, polinomialApproximationπ, propensityπ)

    mteAll = marginalEffect(β10, α, XX, polinomialApproximation, propensity)

    # Pack stuff...
    πSet = Dict(:propensity => propensityπ, 
                :Y => Yπ, 
                :XX => XXπ, 
                :D => Dπ)

    allSet = Dict(:propensity => propensity, 
                :Y => Y, 
                :XX => XX, 
                :D => D)

    return mte, mteAll, πSet, allSet

end



function get_average_treatmet(mte, D, propensity)

    ω_ate = 1
    
    ω_att = (1 .- propensity)./mean(D)

    ω_atu = propensity./(1 - mean(D))

    ATE = mean(mte.*ω_ate)

    ATT = mean(mte.*ω_att)

    ATU = mean(mte.*ω_atu)

    return ATE, ATT, ATU, ω_att, ω_atu

end


function bootstrap_se_model(df, depvar, covariate_names, support)
    
    ATE_list = []; ATU_list = []; ATT_list = []
    for b in 1:1000
        index_b = rand(1:size(df,1),4000)
        df_b = df[index_b,:]
        try
            mte, mteall, πSet, allSet  = get_marginal_response(df_b, depvar, covariate_names, K, support[1], support[2])
            ATE, ATU, ATT = get_average_treatmet(mteall, allSet[:D], allSet[:propensity])
            append!(ATE_list, ATE)
            append!(ATU_list, ATU)
            append!(ATT_list, ATT)            
        catch
            println("Bootstrap error")
        end
    end

    ate_se = std(ATE_list)/sqrt(length(ATE_list)); 
    atu_se = std(ATU_list)/sqrt(length(ATE_list)); 
    att_se = std(ATT_list)/sqrt(length(ATE_list))
    
    return ate_se, atu_se, att_se

end



############################
# Implementation:
############################

# F) With covariates:
#--------------------

df = load_dataset()
model_variables = define_variables()

covariates=[]
append!(covariates,model_variables[:household])
append!(covariates,model_variables[:comunity])
covariate_names = vcat(covariates, [:constant])

covariate_names = vcat(covariates, [:constant])
supp  = [0.1, 1]

#Drop missings:
df = df[:,vcat([model_variables[:depvar_c3], model_variables[:depvar_d2]], 
                            [model_variables[:treatment]], 
                            covariate_names,
                            model_variables[:instruments])]

df = df[all.(!ismissing, eachrow(df)), :]

K = 1

# Hours worked
depvar = model_variables[:depvar_c3]
mte, mteall, πSet, allSet = get_marginal_response(df, depvar, covariate_names, K, supp[1],supp[2]) #ATE, ATU, ATT = get_average_treatmet(mte, πSet[:D], πSet[:propensity])
ATE, ATU, ATT = get_average_treatmet(mteall, allSet[:D], allSet[:propensity])
ate_se, atu_se, att_se = bootstrap_se_model(df, depvar, covariate_names, supp)



# Life Satisfaction
depvar = model_variables[:depvar_d2]
mte, mteall, πSet, allSet = get_marginal_response(df, depvar, covariate_names, K, supp[1],supp[2]) #ATE, ATU, ATT = get_average_treatmet(mte, πSet[:D], πSet[:propensity])
ATE, ATU, ATT = get_average_treatmet(mteall, allSet[:D], allSet[:propensity])
ate_se, atu_se, att_se = bootstrap_se_model(df, depvar, covariate_names, supp)









