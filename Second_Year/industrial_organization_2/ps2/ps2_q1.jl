using Pkg

Pkg.instantiate()


using Distributions
using LinearAlgebra
using DataFrames
using StatFiles
using Tables
using CSV
using Optim
using Random
using Plots
using RDatasets
using StatsPlots


df = DataFrame(CSV.File("ps2_ex1.csv"));

# params...
β = 3
ϕ = 4
δ = .6

# data...
x = df.x
n = df.n

# Generate profits data
π_m = β .* x .- ϕ .- δ .* log.(n) .+ rand(Normal(0,1), size(x,1))


function objective(params, π_m, x, n)

    β = params[1]
    ϕ = params[2]
    δ = params[3]
    resid = π_m .-  β .* x .+ ϕ .+ δ .* log.(n) 
    obj = sum(log.(pdf.(Normal(0,1), resid)))

    return -obj

end

# Anonymous function:
func_anon(p) = objective(p, π_m, x, n)


# Initial parameters:
params_init = rand(3)

# Optimize:
result = optimize(func_anon, params_init, GradientDescent(), Optim.Options(outer_iterations = 1500,
                    iterations= 10000,
                    show_trace=true,
                    show_every=100))

# Obtain Minimizer
param_hat = Optim.minimizer(result)

scatter([3,4,.6],param_hat,  
            label="Parameters", 
            xlabel = "True Parameter", 
            ylabel = "Estimated Parameter",
            legend=:topleft)
plot!(0:0.1:4, 0:0.1:4, label = "45 degree line")

savefig("q1_estimate_parameters.pdf")



