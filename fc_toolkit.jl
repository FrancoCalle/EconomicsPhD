using Pkg

#Install ...
# Pkg.activate(".") 
# Pkg.add("Distributions")
# Pkg.add("StatsBase")
# Pkg.add(["DataFrames","DataFramesMeta","Chain"])
# Pkg.add("Plots")
# Pkg.add("CategoricalArrays")
# Pkg.add("StatFiles")
# Pkg.add("Tables")
# Pkg.add("CSV")
# Pkg.add("Optim")
# Pkg.add("Missings")
# Pkg.instantiate()

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
using Missings

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
    flag_cl::Bool # clusters
    cl::Array # clusters

    # Define constructor function
    function olsRegression(y, x, fe = nothing, cl = nothing, constant = true)
        """
        Inputs: 
        Y: An array of N × 1 dimension (outputs).
        X: An array of N × K dimension (inputs).
        cl: An array of N × missing dimension (clusters)
        """
        
        if isnothing(fe)
            if constant == true
                x = hcat(x,ones(size(x)[1],1))
            else
                x = x
            end
        else
            x = hcat(x, fe)
        end

        xᵀx = transpose(x)*x
        xᵀy = transpose(x)*y
        β = inv(xᵀx)*xᵀy #compute parameters

        isnothing(cl) ? new(β, x, y, false, [cl]) : new(β, x, y, true, cl)
    
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


@doc """
    tsls_regression(y, d, z, x=nothing, fe=nothing, intercept=true)\n
    Inputs
    y: Outcome array N x 1 features
    d: Treatment variable array of N x 1
    z: Instrument vector Array of N x F
    x: Control variables N x K, it can be ::nothing
    fe: Dichotomous Array N x nFixedEffects
    intercept: Bool (default = true)

    Output
    β: Second stage coefficients
    Π: First stage coefficients
""" ->
struct tsls_regression
    
    β::Array{Float64} # coefficient
    x::Array{Float64} # features
    y::Array{Float64} # response
    d::Array{Float64}
    cl::Array # clusters

    function tsls_regression(y, d, z, x=nothing, fe=nothing, cl=nothing, intercept=true)

        n = size(z)[1]
        
        if ~isnothing(x)
            z = hcat(z,x)
            d = hcat(d,x)
        end

        if ~isnothing(fe)
            intercept = false
            z = hcat(z,fe)
            d = hcat(d,fe)
        end

        if intercept == true
            z = hcat(z,ones(n,1))
            d = hcat(d,ones(n,1))
        end
        
        Q_Z, R_Z = qr(z)
        
        ZZ_inv = inv(cholesky(R_Z' * R_Z))
        
        P_Z = z * ZZ_inv * z'
        
        β = inv(d' * P_Z * d) * d' * P_Z * y
        
        Π = z\d
        
        return isnothing(cl) ? new(β, y, d, z, [NaN]) : new(β, y, d, z, cl)

    end
end


function predict_outcome(fit::olsRegression, data = nothing)
    
    isnothing(data) ? fitted = fit.x * fit.β : fitted = data * fit.β

    return(fitted)

end


function predict_outcome(fit::tsls_regression, data = nothing)
    
    isnothing(data) ? fitted = fit.d * fit.β : fitted = data * fit.β

    return(fitted)

end

function predict_outcome(fit::probitModel, data = nothing)
    
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
    
    pscore = predict_outcome(probit)
    
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

    res = y - predict_outcome(fit)
    xᵀx = x'*x

    clusters = unique(cl)
    C = length(clusters)
    R = (C / (C-1)) * ((n-1)/(n-k))

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
    
    Meat = reduce(+, meat_cluster.*R)
    Bread = inv(xᵀx)

    varcov =   Bread' * Meat * Bread
    
    function mat_posdef_fix(X::Matrix; tol = 1e-10) # Thanks Ed.
        if any(diag(X) .< tol)
            e_vals, e_vecs = eigen(Symmetric(X))
            e_vals[e_vals .<= tol] .= tol
            X = e_vecs * Diagonal(e_vals) * e_vecs'
        end
        return X
    end

    vcov_matrix = mat_posdef_fix(varcov)
    
    se = sqrt.(diag(varcov))

    return se

end


function se_homoskedastic(fit::olsRegression)

    x = fit.x; y = fit.y; β = fit.β
    
    N = length(y); K = size(x, 2);

    u = y - predict_outcome(fit) # residuals

    XX_inv = inv(x' * x)

    covar = sum(u.^2) * XX_inv

    covar = covar .* (1 / (N - K)) # dof adjustment

    # Get standard errors, t-statistics, and p-values

    se = sqrt.(covar[diagind(covar)])

    return se

end


function inference(fit::olsRegression)

    # Obtain data parameters
    x = fit.x
    y = fit.y
    β = fit.β
    
    N = length(y); 
    K = size(x, 2);

    u = y - predict_outcome(fit) # residuals

    # Calculate the covariance under homoskedasticity
    if fit.flag_cl == true  # If cluster is passed on struct
        # println("Errors: Clustered")
        se = se_cluster(fit)
    else # If cluster is not pased on struct
        # println("Errors: Homoskedastic")
        se = se_homoskedastic(fit)
    end

    #Get standard errors, t-statistics, and p-values
    # se = sqrt.(covar[diagind(covar)])
    t_stat = β ./ se
    p_val = 2 .* cdf.(TDist(N-K), - abs.(t_stat))
    r2 = 1 - sum(u.^2)/sum((y.-mean(y)).^2)
    # Organize and return output
    output = (β = β, se = se, t = t_stat, p = p_val, r = r2)

    return output

end


@doc """
    select_variables(df::DataFrame, y_name::list, x_name::list, z_name::list)\n
    Inputs \n
    df: DataFrame
    y_name: List of symbols with outcome names. 
    x_name: List of symbols with covariates names. 
    d_name: List of symbols with endogenous. 
    z_name: List of symbols with exogenous. 

    Output \n
    y: N x 1 Array with missings disallowed
    x: N x d_H Array with missings disallowed
    d: N x d_F Array with missings disallowed
""" ->
function select_variables(df::DataFrame, y_name, x_name=nothing, d_name=nothing, z_name=nothing)

    # Join all names in all_names
    all_names = []
    try append!(all_names, y_name) catch end
    try append!(all_names, x_name) catch end
    try append!(all_names, d_name) catch end
    try append!(all_names, z_name) catch end

    for n in 1:length(all_names)

        df = filter(all_names[n] => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)

    end

    y = disallowmissing(Matrix(df[:, y_name]))

    if ~isnothing(d_name)
        d = disallowmissing(Matrix(df[:, d_name]))
    else
        d=nothing
    end

    if ~isnothing(x_name)
        x = disallowmissing(Matrix(Matrix(df[:, x_name])))
    else
        x=nothing
    end

    if ~isnothing(z_name)
        z = disallowmissing(Matrix(df[:, z_name]))
    else
        z=nothing
    end
    
    return y, x, d, z

end
