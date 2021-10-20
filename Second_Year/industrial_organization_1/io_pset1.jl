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
Pkg.instantiate()

#Load packages ...

using CategoricalArrays
using Distributions
using LinearAlgebra
using DataFrames
using StatFiles
# using Chain
using Tables
using CSV
using Optim
using Random

include(joinpath("..", "..","fc_toolkit.jl"))

# Question 2

function labor_demand(k_it, a_it, w_it, p, α)
    
    l = 1/(1-α) * (a_it + log(α) - w_it) + k_it

    return l

end

function production_function(α, ω_ilag)
    
    # Define parameters:
    p = 1 
    ρ = 0.8
    f_γ = Normal(0, 0.5)
    f_ϵ = Normal(0, 0.1)
    f_k = Normal(0, 0.1)
    f_w = Normal(0, 0.5)

    γ_i = rand(f_γ,1)[1]
    ϵ_it = rand(f_ϵ,1)[1]
    k_it = rand(f_k,1)[1]
    w_it = rand(f_w,1)[1]

    ω_it = ρ * ω_ilag + ϵ_it
    a_it = γ_i + ω_it

    l_it = labor_demand(k_it, a_it, w_it, p, α)

    y_it = a_it + α * l_it + (1-α) * k_it

    return y_it, l_it, k_it, ω_it, w_it
end


# B. Construct output, labor and k:

nFirms = 100
T = 50

α = 0.7

y_placeholder = zeros(T*nFirms,1)
l_placeholder = zeros(T*nFirms,1)
k_placeholder = zeros(T*nFirms,1)
ω_placeholder = zeros(T*nFirms,1)
w_placeholder = zeros(T*nFirms,1)
firm_fe = zeros(T*nFirms,1)

for ff in 1:nFirms
    for tt in 1:T
        if tt ==1
            y_placeholder[T*(ff-1) + tt], l_placeholder[T*(ff-1) + tt], k_placeholder[T*(ff-1) + tt], ω_placeholder[T*(ff-1) + tt], w_placeholder[T*(ff-1) + tt] = production_function(α,0)
        elseif tt !=1
            y_placeholder[T*(ff-1) + tt], l_placeholder[T*(ff-1) + tt], k_placeholder[T*(ff-1) + tt], ω_placeholder[T*(ff-1) + tt], w_placeholder[T*(ff-1) + tt] = production_function(α,ω_placeholder[tt-1])
        end
        firm_fe[T*(ff-1) + tt] = ff
    end
end



# C. Now estimate production function using OLS:
X = hcat(l_placeholder,k_placeholder); Y = y_placeholder;
fitOLS = olsRegression(X, Y, nothing, nothing);
res = inference(fitOLS);

println("Addition of participation parameters: ", round(res.β[1] + res.β[2], digits=3), "\n", "β coefficients: ", round.(res.β, digits = 3), " p-values: ", round.(res.p, digits = 3))


# D. Now estimate production function using fixed effects:
# Fixed effects.
FixedEffects = reduce(hcat, [firm_fe .== fe for fe in unique(firm_fe)]);
fitOLS = olsRegression(X, Y, FixedEffects, nothing);
res = inference(fitOLS);

println("Addition of participation parameters: ", round(res.β[1] + res.β[2], digits=3), "\n", "β coefficients: ", round.(res.β[1:2], digits = 3), " p-values: ", round.(res.p[1:2], digits = 3))

# Random effects.

Pkg.add("MixedModels")
using MixedModels

df = DataFrame(:y => y_placeholder[:], :l => l_placeholder[:], :k => k_placeholder[:], :w => w_placeholder[:], :id => CategoricalArray(Int.(firm_fe[:])))

Pkg.add("StatsModels")
using StatsModels
# fm = @formula(y ~ 1 + l + k + (1|firm_fe))
fm1 = fit(MixedModel, @formula(y ~ 1 + l + k + (1|id)), df);
# Basically the same
fm1.β

# E. Is there anything to use as instrument? ...
Pkg.add("Econometrics")
using Econometrics

model = fit(Econometrics,@formula(y ~ k + (l ~ w)),df)

# Found point estimate:
β = tsls_regression(Y, X, w_placeholder, FixedEffects[:,1:end-1])



