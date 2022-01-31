using Pkg

# #Install ...
# Pkg.add("Distributions")
# Pkg.add("StatsBase")
# Pkg.add("Plots")
# Pkg.add("StatFiles")
# Pkg.add("Tables")
# Pkg.add("CSV")
# Pkg.add("Optim")
# Pkg.instantiate()

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
using NLsolve
using ForwardDiff


# Excercise 2

function unpackVariables(dataset)
    
    s_jt = dataset.Shares
    p = dataset.Prices
    x = dataset.x
    z = dataset.z
    m = dataset.market
    j = dataset.Product
    M = maximum(dataset.market);
    J = maximum(dataset.Product);

    return s_jt , p, x , z, m, j, M, J

end


function blpShares(δ_jt, m)

    M = maximum(m)

    s_blp_jt = zeros(size(m))
    
    for ii = 1:M
        s_blp_jt[m .== ii] = exp.(δ_jt[m .== ii])./(1 + sum(exp.(δ_jt[m .== ii])))
    end

    return s_blp_jt
end


function tsls_regression(y, p, z, x=nothing, intercept=true)

    n = size(z)[1]
    
    z = hcat(z,x)
    d = hcat(p,x)

    if intercept==true
        z = hcat(z,ones(n,1))
        d = hcat(d,ones(n,1))
    end

    
    _, R_Z = qr(z)
    
    ZZ_inv = inv(cholesky(R_Z' * R_Z))
    
    P_Z = z * ZZ_inv * z'
    
    β = inv(d' * P_Z * d) * d' * P_Z * y
    
    Π = z\d
    
    return β, Π, d

end


function get_elasticities(α, p, s_jt)

    # Δ = reshape(δ_jt, J, M)
    S = reshape(s_jt, J, M)
    P = reshape(p, J, M)

    η_jkm = zeros(J,J,M)

    for mm = 1:M
        for jj = 1:J
            for kk = 1:J
                if jj == kk
                    η_jkm[jj,kk,mm] = -α .* P[kk,mm] .* (1 - S[kk,mm])
                else
                    η_jkm[jj,kk,mm] = α .* P[kk,mm] .* S[kk,mm]
                end
            end
        end
    end

    return η_jkm
end



function get_marginal_cost(p, η_jkm)

    P = reshape(p, J, M)

    I_J= Matrix(I,J,J)
    mc_jm = zeros(J,M)
    for mm = 1:M 
        mc_jm[:,mm] = P[:,mm] .+ inv(I_J.*η_jkm[:,:,mm]) * P[:,mm]
    end

    return mc_jm

end


# Execute code:...

dataset = DataFrame(CSV.File("ps1_ex3.csv"));

s_jt , p, x , z, m, j, M, J = unpackVariables(dataset);


# Generate Random Mean Utilities ...

δ_jt_lag = rand(M*J);

eps = 1 # tolerance

# Compute contraction ...

while eps > 10^-25

    δ_jt = δ_jt_lag + log.(s_jt) - log.(blpShares(δ_jt_lag,m))
    
    eps = sum((δ_jt - δ_jt_lag).^2)
    
    print("Euclidean Distance: ", eps, "\n")
    
    δ_jt_lag = δ_jt

end

δ_jt = δ_jt_lag


# Just recovered δ's:
histogram(δ_jt)

# Use moment condition and instrumentalize prices:
β, Π, d = tsls_regression(δ_jt, p, z, x, true);
α = -β[1]

ξ = δ_jt - d*β #Estimate unobserved ξ:


# Now compute elasticities:

η_jkm = get_elasticities(α, p, s_jt);

η_jk = zeros(J,J)
for mm = 1:M η_jk += η_jkm[:,:,mm]./M end # Average across markets:

display(η_jk)

CSV.write("Q3_JJ_avg_elasticities.csv", Tables.table(round.(η_jk,digits=3)), writeheader=false) 


# Marginal Cost:

mc_jm = get_marginal_cost(p, η_jkm);
mean(mc_jm, dims=2)

scatter(p, mc_jm[:], xlabel ="Price", ylabel ="Marginal Cost")
savefig("Q3_mc_vs_prices.pdf")


# Part IV of question 2. this one can be tricky ...

# Product j = 1 exits the market, if so, shares change, but by 2 mechanisms
# lower prices, and just higher demand for less products.

# Recall the profit maximization FOC:

# Define the system of equations:



function eq_system(F, x; XX, mc, ξ, β, J = 6)

    # Find Delta Functions:
    F[1:J] = x[J+1:J+J]  .- (x[1:J] .* β[1]  .+ XX .* β[2] .+ 1 .* β[3]  .+ ξ);

    # Find Share Functions:
    F[J+1:J+J] = x[J+J+1:J+J+J] .- exp.(x[J+1:J+J])./(1 + sum(exp.(x[J+1:J+J])))

    # Find Price Functions:
    F[J+J+1:J+J+J] = x[1:J] .- x[1:J]./abs.(β[1].*x[1:J].*(1 .-x[J+J+1:J+J+J])) .- mc

end


# Recursive form wouldn't work, I think we might want to optimize using objective function

ξ = reshape(ξ, J,M)

XX = reshape(d[:,2], J, M)

p_new = zeros(J,M)

shares_new = zeros(J,M)

for mm in 1:M
    anon_eq_system(F, θ) = eq_system(F, θ; 
                                        XX=XX[:,mm], 
                                        mc=mc_jm[:,mm], 
                                        ξ= ξ[:,mm], 
                                        β=β,
                                        J = 6)

    initial_x = rand(18)./100
    initial_F = similar(initial_x)
    df = OnceDifferentiable(anon_eq_system, initial_x, initial_F)
    result = nlsolve(df, initial_x)
    p_new[:,mm] =  result.zero[1:6]
    shares_new[:,mm] =  result.zero[13:18]

end


scatter(p, p_new[:], y_lims=(1.5,5), x_lims=(1.5,5), xlabel ="Original Price", ylabel ="Computed Price",label="Prices")
savefig("Q4_verify_price_optimization.pdf")

scatter(s_jt, shares_new[:], y_lims=(-.1,.5), x_lims=(-.1,.5), xlabel ="Original Share", ylabel ="Computed Share",label="Shares")
savefig("Q4_verify_shares_optimization.pdf")



# Now let's do simulation ...

p_new = zeros(5,M)
shares_new = zeros(5,M)

for mm in 1:M
    anon_eq_system(F, θ) = eq_system(F, θ; 
                                        XX=XX[2:end,mm], 
                                        mc=mc_jm[2:end,mm], 
                                        ξ= ξ[2:end,mm], 
                                        β=β,
                                        J = 5)

    initial_x = rand(15)./100
    initial_F = similar(initial_x)
    df = OnceDifferentiable(anon_eq_system, initial_x, initial_F)
    result = nlsolve(df, initial_x)
    p_new[:,mm] =  result.zero[1:5]
    shares_new[:,mm] =  result.zero[11:15]

end


mean(reshape(p, J,M ), dims=2)

# Average prices and shares after simulation
print("Average Prices: ", round.(mean(p_new, dims=2), digits= 5))
print("Average Shares: ",round.(mean(shares_new, dims=2), digits= 5))


# Part 5:

function profit_maximization(p, mc, s)

    π = (p .- mc).*s

    return π

end

function welfare(p, XX, ξ, s, β)

    utility = β[1].*p .+ β[2].*XX .+ β[3] .+ ξ
    utility = utility .* s
    total_welfare = sum(utility .* s, dims = 1)

    return total_welfare

end

profit_0 = profit_maximization(reshape(p,J,M), mc_jm, reshape(s_jt,J,M))

profit_1 = profit_maximization(p_new, mc_jm[2:end,:], shares_new)

increase_profits = (profit_1 .- profit_0[2:end,:])./profit_0[2:end,:]


print("Average Increase in Profits: ", round.(mean(increase_profits, dims=2), digits=5))


# Consumer welfare:

total_welfare0 = welfare(reshape(p, J,M), XX, ξ, reshape(s_jt, J,M), β)

total_welfare1 = welfare(p_new, XX[2:end,:], ξ[2:end,:], shares_new, β)

R = abs(min(minimum(total_welfare1), minimum(total_welfare0)))

increase_welfare = ((total_welfare1.+R) .- (total_welfare0.+R))./(total_welfare0.+R)

print("Average Welfare Increase: ", mean(increase_welfare, dims=2))

histogram(increase_welfare[:], x_lims=(-1,1))


