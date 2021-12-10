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
using MyMethods

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



# Functions for part H & I:


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
                    sense="Max", decreasing=false, prior_bound=(1,-9))

    # Data parameters
    K = length(gamma_1IV) - 1

    # initialize model
    m = Model(GLPK.Optimizer)

    # initialize variables
    @variable(m, prior_bound[1] >= theta[1:(K+1), 1:2] >=  prior_bound[2]) # bounded by 0 and 1

    # set constraints
    @constraint(m, IV, sum(theta[:,1].*gamma_1IV') +
        sum(theta[:,2].*gamma_0IV')  == b_IV)

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


df = load_dataset()
model_variables = define_variables()
covariate_names = [:constant]


# Part I: Parameterize the MTR function using the Bernstein polynomials, omiting
# covariates, and re-compute bounds on the ATU using the same prior
# outcome bounds as in (h).


# Get Bernstein poly basis terms and compute the gamma
b_IV = -3.52
b_TSLS = -3.5
K = 30
nMC = size(df,1)
u = collect(1:nMC)./(nMC+1) #rand(Uniform(),nMC)
w1_IV, w0_IV, w1_TSLS, w0_TSLS, w1_ATT, w1_ATU = get_weights(df)

# We have 3: ATU 
p_bounds = zeros(2, 2, K)
pb = (0,-9)
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
                sense = "Min", decreasing = dec, prior_bound=pb)

        # upper bound
        p_bounds[dec+1, 2, k] = get_bound(b_IV, b_TSLS,
                gamma_1IV, gamma_0IV,
                gamma_1TSLS, gamma_0TSLS,
                gamma_1ATU, gamma_0ATU,
                sense = "Max", decreasing = dec, prior_bound=pb)
    end
end


pyplot(size=(500,300), leg=true);
_x = collect(1:K)

# parametric bounds
plot(_x, p_bounds[1,1,:], line = (:line, :line, 0.8, 1, :blue), xlabel = "Polynomial degree", ylabel = "Upper and Lower Bounds", label=nothing)
plot!(_x, p_bounds[1,1,:], seriestype = :scatter, markershape=:circle, markersize=4, color=:blue, label="ATU Bounds")
plot!(_x, p_bounds[1,2,:], line = (:line, :line, 0.8, 1, :blue), label=nothing)
plot!(_x, p_bounds[1,2,:], seriestype = :scatter, markershape=:circle, markersize=4, color=:blue, label=nothing)


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
            sense = "Min", decreasing = dec, prior_bound=pb)
    
    # upper bound
    np_bounds[dec+1, 2] = get_bound(b_IV, b_TSLS,
            gamma_1IV, gamma_0IV,
            gamma_1TSLS, gamma_0TSLS,
            gamma_1ATU, gamma_0ATU,
            sense = "Max", decreasing = dec, prior_bound=pb)
end

# cast np_bounds to dimension of p_bounds
np_bounds = reshape(repeat(np_bounds, K)', (2, 2, K));

plot!(_x, np_bounds[1,1,:], line = (:line, :dot, 0.8, 1, :red), label="Nonparametric")
plot!(_x, np_bounds[2,1,:], line = (:line, :dot, 0.8, 1, :red), label="")
savefig("Q1_PH_ATU_noparam_bounds_work_hours.pdf")







# Life Satisfaction:
#--------------------

# Get Bernstein poly basis terms and compute the gamma
b_IV = 1.95
# b_TSLS = -3.5
K = 30
nMC = size(df,1)
u = rand(Uniform(),nMC)
w1_IV, w0_IV, w1_TSLS, w0_TSLS, w1_ATT, w1_ATU = get_weights(df)

# We have 3: ATU 
p_bounds = zeros(2, 2, K)
pb = (3,0)
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
                sense = "Min", decreasing = dec, prior_bound=pb)

        # upper bound
        p_bounds[dec+1, 2, k] = get_bound(b_IV, b_TSLS,
                gamma_1IV, gamma_0IV,
                gamma_1TSLS, gamma_0TSLS,
                gamma_1ATU, gamma_0ATU,
                sense = "Max", decreasing = dec, prior_bound=pb)
    end
end


pyplot(size=(500,300), leg=true);
_x = collect(1:K)

# parametric bounds
plot(_x, p_bounds[1,1,:], line = (:line, :line, 0.8, 1, :blue), xlabel = "Polynomial degree", ylabel = "Upper and Lower Bounds", label=nothing)
plot!(_x, p_bounds[1,1,:], seriestype = :scatter, markershape=:circle, markersize=4, color=:blue, label="ATU Bounds")
plot!(_x, p_bounds[1,2,:], line = (:line, :line, 0.8, 1, :blue), label=nothing)
plot!(_x, p_bounds[1,2,:], seriestype = :scatter, markershape=:circle, markersize=4, color=:blue, label=nothing)


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
            sense = "Min", decreasing = dec, prior_bound=pb)
    
    # upper bound
    np_bounds[dec+1, 2] = get_bound(b_IV, b_TSLS,
            gamma_1IV, gamma_0IV,
            gamma_1TSLS, gamma_0TSLS,
            gamma_1ATU, gamma_0ATU,
            sense = "Max", decreasing = dec, prior_bound=pb)
end

# cast np_bounds to dimension of p_bounds
np_bounds = reshape(repeat(np_bounds, K)', (2, 2, K));

plot!(_x, np_bounds[1,1,:], line = (:line, :dot, 0.8, 1, :red), label="Nonparametric")
plot!(_x, np_bounds[2,1,:], line = (:line, :dot, 0.8, 1, :red), label="")
savefig("Q1_PH_ATU_noparam_bounds_life_satisfaction.pdf")


