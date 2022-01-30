# Problem Set 1; Chuhan and Calle:
#---------------------------------

using Pkg

#Install ...
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add("Plots")
Pkg.add("StatFiles")
Pkg.add("Tables")
Pkg.add("CSV")
Pkg.add("Optim")
Pkg.instantiate()

#Load packages ...

using Distributions
using LinearAlgebra
using DataFrames
using StatFiles
using Chain
using Tables
using CSV
using Optim
using Random
using Plots


# Excercise 2

struct ModelData
    Y::Array
    X::Array
    XX::Array
    D::Array
    N::Int
    J::Int
end

struct Parameters
    β::Array
    Γ::Array
end

function unpackVariables(dataset)
    
    Y = dataset.choice
    X = zeros(size(dataset,1), 3)
    D = zeros(size(dataset,1), 2)

    X[:,1] = dataset[:,Symbol("x.1")]
    X[:,2] = dataset[:,Symbol("x.2")]
    X[:,3] = dataset[:,Symbol("x.3")]

    # Rearange x's
    J = maximum(Y)
    XX = zeros(3,J)
    for ii = 1:J XX[:,ii] = X[Y .== ii,:][1,:] end
    
    D[:,1] = dataset[:,Symbol("d.1")]
    D[:,2] = dataset[:,Symbol("d.2")]

    return Y , X, XX , D
end


# First work out the model with fake data...

function fakeUtilityFunction(md::ModelData, mp::Parameters)

    XX = md.XX; D = md.D
    N = md.N; J = md.J

    β = mp.β; Γ = mp.Γ

    # Generate Gumbel dist
    ϵ = rand(Gumbel(0,1), N, J)
    
    # Say ξ ∼ N(0,.2)
    ξ = rand(Normal(0,.2), J)
    ξ[end] = 0

    # Generate utilities.
    u_ij = zeros(N, J)
    Y = zeros(N,1)

    # Compute Mean Utility:
    δ = XX'*β  .+  ξ

    for ii in 1:N 
        u_ij[ii,:] =  δ .+ (D[ii,:]' * Γ * XX)' .+ ϵ[ii,:]
        Y[ii] = argmax(u_ij[ii,:])
    end

    return Y, δ

end


function UtilityFunction(parameters, md::ModelData)
    
    Y = md.Y
    XX = md.XX; D = md.D; N = md.N; J = md.J

    Γ = reshape(parameters[1:6], 2,3)
    δ = zeros(J)
    δ[1:J-1] = parameters[7:end]

    # Generate utilities.
    u_ij = zeros(N, J)
    pr_ij = zeros(N, J)

    for ii in 1:N 
        u_ij[ii,:] =  δ .+ (D[ii,:]' * Γ * XX)' 
    end

    pr_ij = exp.(u_ij)
    pr_ij = pr_ij./sum(pr_ij, dims=2)
    pr_i = Array([pr_ij[ii, Y[ii]] for ii = 1:N])

    return pr_i

end

function objectiveFunction(parameters, md::ModelData)

    pr_i = UtilityFunction(parameters, md)
    logLik = sum(log.(pr_i))
    
    return -logLik

end



# Run stuff: 

dataset = DataFrame(CSV.File("ps1_ex2.csv"))

Y, X, XX, D = unpackVariables(dataset)

Y_set = sort(unique(dataset.choice))

J = maximum(Y_set)

N = size(Y,1)

md = ModelData(Y, X, XX, D, N, J)

mp = Parameters(rand(3), rand(2,3))


# Now Generate fake data based on previous information:

Y_new, δ = fakeUtilityFunction(md, mp)

md_new = ModelData(Int.(Y_new), X, XX, D, N, J)


# Estimate θ = (β, Γ)
param_init = rand(J+6-1) # -1 because of outside option

obj_function(x) = objectiveFunction(x, md_new)

result = optimize(obj_function, param_init, LBFGS())
        
params_hat = Optim.minimizer(result)

Γ_hat = reshape(params_hat[1:6],2,3)

δ_hat = params_hat[7:end]

scatter(Γ_hat[:], mp.Γ[:], xlims = (.3,.8), ylims = (.3,.8), label="parameters: Γ")
plot!(0:0.1:1,0:0.1:1, label="45 degree line")
savefig("Q2_gamma_fit.pdf")

scatter(δ, δ_hat)
plot!(minimum(δ):0.1:maximum(δ),minimum(δ):0.1:maximum(δ), label="45 degree line")
savefig("Q2_mean_utility_fit.pdf")

# We reconver parameters!! 



# Now try with true picks...

# Estimate θ = (β, Γ)
param_init = rand(J+6-1) # -1 because we eliminate mean utility of outside option

obj_function(x) = objectiveFunction(x, md)

result = optimize(obj_function, param_init, LBFGS())
        
params_hat = Optim.minimizer(result)

Γ_hat = reshape(params_hat[1:6],2,3)

δ_hat = params_hat[7:end]

histogram(δ_hat, label="Mean Utility: δⱼ")
savefig("Q2_mean_utility_distribution.pdf")

# bar(Γ_hat[1,:], label="D1 vs X")
# bar(Γ_hat[2,:], label="D2 vs X")


# Use moment condition E[X ξ] = 0

β_hat = XX'\δ
print("The parameters β are: ",β_hat)
