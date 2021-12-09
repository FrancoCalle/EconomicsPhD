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


# Part E: Estimate ATE, ATU, and ATT using MTE.
# - Restrict attention to parametric specifications of MTR
# - 1st estimating the treatment selection equation in (2) as a probit model to obtain estimates of the propensity score
# - 2nd modeling 𝐾(𝑝) as a polynomial in 𝑝 of degree k and estimating the outcome equation

# This one uses a probit built in function, planning to change it later: 

function get_propensity_score(df,covs)

    if length(covs) > 1
        reg = @formula(connected ~ treat_full + treat_low + treat_med + female + base_age + educ + bank + housing + asset_value + energyspending + busia + market + transearly + connected_rate + population)
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

# Part E: Estimate ATE, ATU, and ATT using MTE.
# - Restrict attention to parametric specifications of MTR
# - 1st estimating the treatment selection equation in (2) as a probit model to obtain estimates of the propensity score
# - 2nd modeling 𝐾(𝑝) as a polynomial in 𝑝 of degree k and estimating the outcome equation


# F) With covariates:
#--------------------
K = 2
covariate_names = vcat(covariates, [:constant])
support  = (0.01, 1)

# Hours worked
depvar = model_variables[:depvar_c3]
mte, mteall, πSet, allSet = get_marginal_response(df, depvar, covariate_names, K, support[1],support[2]) #ATE, ATU, ATT = get_average_treatmet(mte, πSet[:D], πSet[:propensity])
ATE, ATU, ATT = get_average_treatmet(mteall, allSet[:D], allSet[:propensity])
ate_se, atu_se, att_se = bootstrap_se_model(df, depvar, covariate_names, support)

# Life Satisfaction
depvar = model_variables[:depvar_d2]
mte, mteall, πSet, allSet = get_marginal_response(df, depvar, covariate_names, K, support[1],support[2]) #ATE, ATU, ATT = get_average_treatmet(mte, πSet[:D], πSet[:propensity])
ATE, ATU, ATT = get_average_treatmet(mteall, allSet[:D], allSet[:propensity])
ate_se, atu_se, att_se = bootstrap_se_model(df, depvar, covariate_names, support)


histogram(allSet[:propensity][allSet[:D][:].==1],  bins =20 , fillalpha=0.2)
histogram!(allSet[:propensity][allSet[:D][:].==0], bins =20 ,fillalpha=0.2)

scatter(allSet[:propensity],mteall,fillalpha=0.2)



# G) Now consider policy where effective connection price is changed to $200
# - Construct estimate of per person policy relevant treatment relative to status quo
# - use two outcomes in (a) and parametric, point identified specifications of mtr as in e.

# Small tweak, instead of dichotomous variables I'll transform it in effective connection price:
# - MTE will remain the same, we will just change propensity score and compute PRTE

function compute_prte(df, depvar, covariate_names, effective_cost)
    
    K = 2

    support  = (0.01, 1)

    _, mteall, _, _ = get_marginal_response(df, depvar, covariate_names, K, support[1], support[2]) #ATE, ATU, ATT = get_average_treatmet(mte, πSet[:D], πSet[:propensity])

    probit = glm(@formula(connected ~ effective_connection_price), df, Binomial(), ProbitLink())

    propensity_0 = predict(probit,df)

    df_new = df; df_new[:,:effective_connection_price] .= effective_cost

    propensity_1 = predict(probit,df_new)

    # Compute mte on identified segment:

    PRTE = mean(mteall.*(propensity_0.-propensity_1)./(mean(propensity_1)-mean(propensity_0)))

    return PRTE

end


function bootstrap_prte(df, depvar, covariate_names, effective_cost)
    
    PRTE_list = [];
    for b in 1:1000
        index_b = rand(1:size(df,1),4000)
        df_b = df[index_b,:]
        try
            PRTE = compute_prte(df_b, depvar, covariate_names, effective_cost)
            append!(PRTE_list, PRTE)            
        catch
            println("Bootstrap error")
        end
    end

    prte_se = std(PRTE_list)/sqrt(length(PRTE_list)) 

    return prte_se

end



# Hours Work
covariate_names = [:constant]
depvar = model_variables[:depvar_c3]
PRTE = compute_prte(load_dataset(), depvar, covariate_names, 200)
PRTE_se = bootstrap_prte(load_dataset(), depvar, covariate_names, 200)

# Life Satisfaction
covariate_names = [:constant]
depvar = model_variables[:depvar_d2]
PRTE = compute_prte(load_dataset(), depvar, covariate_names, 200)
PRTE_se = bootstrap_prte(load_dataset(), depvar, covariate_names, 200)




# Functions for part H & I:

using MyMethods


#From here I will use the notation from the jupyter notebook applied to my problem:
# Z ∈ {low, med, high}
# p 

function get_weights(df)
    
    Z = Array(df[:,model_variables[:instruments]])
    D = Array(df[:,model_variables[:treatment]])
    
    ind_1 = copy(Z)
    for ii in 1:size(ind_1,1)
        for jj in 1:2
            if ind_1[ii,end-jj+1] == 1 
                ind_1[ii,end-jj] = 1
            end
        end
    end

    z = Array(1:3)
    # probability of getting treatment z ∈ {1,2,3}
    p_z = mean(Z,dims=1)[:]
    p_z = p_z./sum(p_z)
    
    # Propensity: p = P[D = 1 | Z = z]
    p = [mean(D[Z[:,ii] .== 1]) for ii = 1:3]

    Cov_DZ = (sum(p_z.*z.*p)) - sum(p.*p_z)*sum(z.*p_z) # Define population moments
    s_IV = ((z .- sum(z.*p_z))./Cov_DZ)' # s(d,z)
    w1_IV = s_IV .* ind_1 # weight for d = 1
    w0_IV = s_IV .* (1 .- ind_1) # weight for d = 0

    # TSLS slope
    ZD = hcat(p_z, p_z.*p)' # Define population moments
    inv_ZZ = I.*3
    FS = inv_ZZ * ZD # first stage
    s_TSLS = (inv(FS * ZD') * ZD * inv_ZZ)[2,:]' # s(d,z)
    w1_TSLS = s_TSLS .* ind_1 # weight for d = 1
    w0_TSLS = s_TSLS .* (1 .- ind_1) # weight for d = 0

    # w0_ATT and w0ATU a la Mogstad and Torgovitsky:
    w1_ATT = ind_1./mean(D)
    w1_ATU  = (1 .- ind_1)/(1-mean(D))

    return w1_IV, w0_IV, w1_TSLS, w0_TSLS, w1_ATT, w1_ATU

end


function get_gamma(B, w1, w0)
    
    # Get γ using some polynomial (Bernstein or Non parametrics)
    
    Bw1 = B .* mean(w1, dims=2)
    
    Bw0 = B .* mean(w0, dims=2)
    
    return mapslices(mean, Bw1, dims=1), mapslices(mean, Bw0, dims=1)

end



function get_bound(b_IV, b_TSLS,
                    gamma_1IV, gamma_0IV,
                    gamma_1TSLS, gamma_0TSLS,
                    gamma_1ATT, gamma_0ATT;
                    sense="Max", decreasing=false)

    # Data parameters
    K = length(gamma_1IV) - 1

    # initialize model
    m = Model(GLPK.Optimizer)

    # initialize variables
    @variable(m, 1 >= theta[1:(K+1), 1:2] >= 0) # bounded by 0 and 1

    # set constraints
    @constraint(m, IV, sum(theta[:,1].*gamma_1IV') +
        sum(theta[:,2].*gamma_0IV')  == b_IV)
    # @constraint(m, TSLS, sum(theta[:,1].*gamma_1TSLS') +
    #     sum(theta[:,2].*gamma_0TSLS') == b_TSLS)

    # Restrict to decreasing MTRs
    if decreasing
        @constraint(m, decreasing[j=1:K,s=1:2], theta[j,s] >= theta[j+1,s])
    end

    # define objective 
    if sense == "Max"
    @objective(m, Max, sum(theta[:,1].*gamma_1ATT') +
            sum(theta[:,2].*gamma_0ATT')) # upper bound
    elseif sense == "Min"
    @objective(m, Min, sum(theta[:,1].*gamma_1ATT') +
        sum(theta[:,2].*gamma_0ATT')) # upper bound
    end

    # solve model
    MOI.set(m, MOI.Silent(), true)
    optimize!(m)
    bound = objective_value(m)

    # Return bound
    return bound
end


# Part I: Parameterize the MTR function using the Bernstein polynomials, omiting
# covariates, and re-compute bounds on the ATU using the same prior
# outcome bounds as in (h).


# Get Bernstein poly basis terms and compute the gamma
b_IV = -2.6
b_TSLS = -3.5
K = 30
nMC = length(mte)
u = rand(Uniform(),nMC)
w1_IV, w0_IV, w1_TSLS, w0_TSLS, w1_ATT, w1_ATU = get_weights(df)

# We have 3: ATU 
p_bounds = zeros(2, 2, K)
for k in 1:K
    # Get Bernstein poly basis terms and compute the gamma
    Bernstein = MyMethods.get_basis(u, "Bernstein", k, nothing)
    gamma_1IV, gamma_0IV = get_gamma(Bernstein, w1_IV, w0_IV)
    gamma_1TSLS, gamma_0TSLS = get_gamma(Bernstein, w1_TSLS, w0_TSLS)
    gamma_1ATU, gamma_0ATU = get_gamma(Bernstein, .-w1_ATU, w1_ATU)

    # Compute the parametric bounds
    for dec in (false,true)
        # lower bound
        p_bounds[dec+1, 1, k] = get_bound(b_IV, b_TSLS,
                gamma_1IV, gamma_0IV,
                gamma_1TSLS, gamma_0TSLS,
                gamma_1ATU, gamma_0ATU,
                sense = "Min", decreasing = dec)

        # upper bound
        p_bounds[dec+1, 2, k] = get_bound(b_IV, b_TSLS,
                gamma_1IV, gamma_0IV,
                gamma_1TSLS, gamma_0TSLS,
                gamma_1ATU, gamma_0ATU,
                sense = "Max", decreasing = dec)
    end
end


pyplot(size=(500,300), leg=true);
_x = collect(1:K)

# parametric bounds
plot(_x, p_bounds[1,1,:], line = (:line, :line, 0.8, 1, :blue), xlabel = "Polynomial degree", ylabel = "Upper and Lower Bounds", label=nothing)
plot!(_x, p_bounds[1,1,:], seriestype = :scatter, markershape=:circle, markersize=4, color=:blue, label="ATU Bounds")
plot!(_x, p_bounds[1,2,:], line = (:line, :line, 0.8, 1, :blue), label=nothing)
plot!(_x, p_bounds[1,2,:], seriestype = :scatter, markershape=:circle, markersize=4, color=:blue, label=nothing)
savefig("Q1_PH_ATU_bounds.pdf")




# Part H: Pick prior lower and upper bounds for the two outcomes in A and  discuss choices.
# - Compute nonparametric bounds on the ATU that use all info in E[Y_i | D_i, Z_i].

p = [0.16216216216216217, 0.319060773480663, 0.9473684210526315]
CSplines = MyMethods.get_basis(u, "CSplines2", 6, p) 
gamma_1IV, gamma_0IV = get_gamma(CSplines, w1_IV, w0_IV)
gamma_1TSLS, gamma_0TSLS = get_gamma(CSplines, w1_TSLS, w0_TSLS)
gamma_1ATU, gamma_0ATU = get_gamma(CSplines, .-w1_ATU, w1_ATU)

# Compute the nonparametric bounds
np_bounds = zeros(2,2)
for dec in (false,true)
    # lower bound
    np_bounds[dec+1, 1] = get_bound(b_IV, b_TSLS,
            gamma_1IV, gamma_0IV,
            gamma_1TSLS, gamma_0TSLS,
            gamma_1ATU, gamma_0ATU,
            sense = "Min", decreasing = dec)
    
    # upper bound
    np_bounds[dec+1, 2] = get_bound(b_IV, b_TSLS,
            gamma_1IV, gamma_0IV,
            gamma_1TSLS, gamma_0TSLS,
            gamma_1ATU, gamma_0ATU,
            sense = "Max", decreasing = dec)
end

# cast np_bounds to dimension of p_bounds
np_bounds = reshape(repeat(np_bounds, K)', (2, 2, K));

plot!(_x, np_bounds[1,1,:], line = (:line, :dot, 0.8, 1, :red), label="Nonparametric")
plot!(_x, np_bounds[2,1,:], line = (:line, :dot, 0.8, 1, :red), label="")

savefig("Q1_PH_ATU_noparam_bounds.pdf")













