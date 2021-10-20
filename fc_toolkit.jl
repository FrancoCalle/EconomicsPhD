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

using Distributions
using LinearAlgebra
# using StatsBase
using DataFrames
# using Plots
# using CategoricalArrays
using StatFiles
using Chain
using Tables
using CSV
using Optim
using Random


@doc """
    Inputs\n
    x: Array of N x K features \n
    y: Outcome Array of N x 1 \n
    fe: Binary array of N x nFixedEffects \n
    cl: Categorical array with cluster classification N x 1  \n

    Output \n
    minKidx: minimum K indexes
""" ->
struct olsRegression
    
    β::Array{Float64} # coefficient
    x::Array{Float64} # features
    y::Array{Float64} # response
    cl::Array # clusters

    # Define constructor function
    function olsRegression(x::Array{Float64}, y::Array{Float64}, fe = nothing, cl = nothing)
        """
        Inputs: 
        X: An array of N × K dimension (inputs).
        Y: An array of N × 1 dimension (outputs).
        cl: An array of N × missing dimension (clusters)
        """
        
        if isnothing(fe)
            x = hcat(x,ones(size(x)[1],1))
        else
            x = hcat(x, fe)
        end

        xᵀx = transpose(x)*x
        xᵀy = transpose(x)*y
        β = inv(xᵀx)*xᵀy #compute parameters

        isnothing(cl) ? new(β, x, y, [NaN]) : new(β, x, y, cl)
    
    end

end 


@doc """
    Inputs\n
    x: Array of N x K features \n
    d: Outcome Array of N x 1 \n

    Output \n
    minKidx: minimum K indexes
""" ->
struct probitModel
    
    θ::Array{Float64} # coefficient
    μ::Float64 # Dist. expected value
    σ::Float64 # Dist. variance
    x::Array{Float64} # features
    d::Array{Float64} # dichotomous outcome

    # Define constructor function
    function probitModel(x, d)
        """
        Inputs: 
        x: An array of N × K dimension (inputs).
        d: An dichotomous array of N × 1 dimension (outputs).
        """
        
        x = hcat(x,ones(size(x)[1],1))
        
        param_init  = zeros(size(x,2) + 2,1) .+ 0.01 

        sigmoid(x) = 1/(1+exp(x))

        function objective(params, x, d)
        
            σ = sigmoid(params[end])
            μ = params[end-1]
            θ = params[1:end-2]
            yhat = x*θ
            Φ = Normal(μ,σ)
            log_p = d.*log.(cdf.(Φ, yhat)) + (1 .- d).*log.(1.0 .- cdf.(Φ, yhat))
            log_L = reduce(+,log_p)
            return -log_L
        
        end        

        f(θ) = objective(θ, x, d)
        
        result = optimize(f, param_init, LBFGS())
        
        params_hat = Optim.minimizer(result)

        new(params_hat[1:end-2], params_hat[end-1], sigmoid(params_hat[end]), x, d)

    end

end 


function tsls_regression(y, d, z, x=nothing)
 
    if ~isnothing(x)
        x = hcat(x,ones(size(x)[1],1))  
        z = hcat(z,x)
        d = hcat(d,x)
    end
    
    zᵀz = transpose(z)*z
    zᵀd = transpose(z)*d
    Π = inv(zᵀz) * zᵀd
    pred = transpose(Π)*transpose(z)
    β = inv(pred*d)*pred*y

    return β

end


function predict(fit::olsRegression, data = nothing)
    
    isnothing(data) ? fitted = fit.x * fit.β : fitted = data * fit.β

    return(fitted)

end

function predict(fit::probitModel, data = nothing)
    
    isnothing(data) ? fitted = cdf.(Normal(fit.μ, fit.σ),fit.x * fit.θ) : fitted = cdf.(Normal(fit.μ, fit.σ), data * fit.θ)

    return(fitted)

end

@doc """
    Inputs\n
    k: Number of neighbors to use\n
    pscorei: Propensity score for individual i, should be float\n
    neighborhood: Array that contains pscore for all individuals we are comparing\n

    Output \n
    minKidx: minimum K indexes
""" ->
function kNearest(k::Int, pscorei::Float64, neighborhood::Array{Float64})

    #Calculate euclidean distance:
    distance = (pscorei .- neighborhood).^2 

    Random.seed!(1234); distance = distance.* 100000 .+ shuffle(1:size(distance)[1]) ./ 10000 ; # Random ties elimination
    
    dict = Dict(distance .=> 1:size(distance)[1])
    
    minK = sort(distance)[1:k]
    
    minKidx = [dict[ii] for ii in minK]

    return minKidx

end



@doc """
    Inputs \n
    x: Covariates used to compute score P[D = 1 | X] \n
    d: Treatment variable \n
    y: Outcome variable \n
    k: Number of neighbors \n

    Output \n
    ATE: Average Treatment Effect \n
    ATT: Average Treatment on the Treated 
""" ->
function propensityScoreMatching(x, y, d, k)

    probit = probitModel(x,d)
    
    pscore = predict(probit)
    
    y_1 = y[vec(d .== 1)]
    
    y_0 = y[vec(d .== 0)]
    
    pscore_1 = pscore[vec(d .== 1)]
    
    pscore_0 = pscore[vec(d .== 0)]

    y_cf = zeros(size(y,1),1)

    for ii in 1:size(pscore,1)
        
        if d[ii] == 1 # If treated, compare to untreateds...
            k_index = kNearest(k, pscore[ii], pscore_0) 
            y_mean = mean(y_0[k_index])
        else
            k_index = kNearest(k, pscore[ii], pscore_1)
            y_mean = mean(y_1[k_index])
        end
        y_cf[ii] = y_mean
    end

    ATE = (sum(y[vec(d .== 1)]) + sum(y_cf[vec(d .== 0)]) - sum(y[vec(d .== 0)]) - sum(y_cf[vec(d .== 1)]))/size(y,1)
    ATT = mean(y[vec(d .== 1)]) - mean(y_cf[vec(d .== 1)])

    return ATE, ATT

end



function se_cluster(fit::olsRegression)

    # Obtain data parameters
    x = fit.x
    y = fit.y
    cl = fit.cl
    
    n = size(y, 1); 
    k = size(x, 2);

    res = y - predict(ols)
    xᵀx = x'*x

    clusters = unique(cl)
    C = length(clusters)

    function clust_residuals(c)

        cindex = findall(cl .== c) # Find index of observations that belong to cluster c 
        Xc = Matrix(x[cindex,:])   # Convert to matrix
        resc = Matrix(res[cindex,:]) 
        meat = Xc'*resc*resc'*Xc

        return meat
    end

    meat_cluster = broadcast(
                            c -> clust_residuals(c), 
                            clusters
                            )
    
    Meat = reduce(+, meat_cluster)/n
    Bread = inv(xᵀx)

    varcov = (C*(n-1) / ((C-1)*(n-k))) * n * Bread' * Meat * Bread
    se_cluster = sqrt.(diag(varcov))

    return se_cluster

end


function inference(fit::olsRegression)

    # Obtain data parameters
    x = fit.x
    y = fit.y
    β = fit.β
    
    N = length(y); 
    K = size(x, 2);

    # Calculate the covariance under homoskedasticity
    # if se_option = 'homoskedacity'
    u = y - predict(fit) # residuals
      xᵀx = inv(x' * x)
      covar = sum(u.^2) * xᵀx
      covar = covar .* (1 / (N - K)) # dof adjustment

    #Get standard errors, t-statistics, and p-values
    se = sqrt.(covar[diagind(covar)])
    # se = se_cluster(fit)
    t_stat = β ./ se
    p_val = 2 * cdf.(Normal(), -abs.(t_stat))
    r2 = 1 - sum(u.^2)/sum((y.-mean(y)).^2)
    # Organize and return output
    output = (β = β, se = se, t = t_stat, p = p_val, r = r2)

    return output

end


function select_variables(df::DataFrame, y_name, x_name, d_name=nothing)

    df = filter(y_name[1] => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)

    for n in 1:length(x_name)

        df = filter(x_name[n] => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)

    end

    y = Matrix(df[:, y_name])

    x = Matrix(df[:, x_name])

    if ~isnothing(d_name)
        d = Matrix(df[:, d_name])
    else
        d = nothing
    end
    
    return y, x, d 

end
