using Pkg

#Install ...
Pkg.activate(".") 
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add(["DataFrames","DataFramesMeta","Chain"])
Pkg.add("Plots")
Pkg.add("CategoricalArrays")
Pkg.add("StatFiles")
Pkg.add("CSV")
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
using CSV


struct olsRegression
    
    β::Array{Float64} # coefficient
    x::Array{Float64} # features
    y::Array{Float64} # response
    cl::Array # clusters

    # Define constructor function
    function olsRegression(x, y, fe = nothing, cl = nothing)
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


function predict(fit::olsRegression, data = nothing)
    
    isnothing(data) ? fitted = fit.x * fit.β : fitted = data * fit.β

    return(fitted)

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
    #   u = y - predict(fit) # residuals
    #   xᵀx = inv(x' * x)
    #   covar = sum(u.^2) * xᵀx
    #   covar = covar .* (1 / (N - K)) # dof adjustment

    # Get standard errors, t-statistics, and p-values
    # se = sqrt.(covar[diagind(covar)])
    se = se_cluster(fit)
    t_stat = β ./ se
    p_val = 2 * cdf.(Normal(), -abs.(t_stat))

    # Organize and return output
    output = (β = β, se = se, t = t_stat, p = p_val)

    return output

end


function select_variables(df::DataFrame, y_name, x_name)

    df = filter(y_name[1] => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)

    for n in 1:length(x_name)

        df = filter(x_name[n] => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df)

    end

    y = Matrix(df[:, y_name])

    x = Matrix(df[:, x_name])
    
    return y, x 

end

forcedcoexistence_webfinal = DataFrame(load("11423_Data_and_Programs/forcedcoexistence_webfinal.dta"));

regressionDataset = subset(forcedcoexistence_webfinal, :year => ByRow(year -> year ==(2000)))

regressionDataset[:,:clust_id] = string.(lpad.(regressionDataset[:,:eaid],3,'0')) .* string.(lpad.(regressionDataset[:,:statenumber],2,'0'))

# Replicate Table 3:
# Define dependent variable:
y_name = [:logpcinc]

# MODEL 1
#---------------------------------------------------------
# Define `independent' variables
x_name = [:FC, :HC]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,nothing,regressionDataset[:, :clust_id])
# Get standard errors:
results1 = inference(ols) # This is assuming homoskedacity...
results1.t

# MODEL 2 [Check]
#---------------------------------------------------------
# Define `independent' variables
x_name = [:FC, :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,nothing,regressionDataset[:, :clust_id])
results2 = inference(ols) # This is assuming homoskedacity...
results2.β
results2.t

# MODEL 3 [Check]
#---------------------------------------------------------
# Define `independent' variables
x_name = [:FC, :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm, :ea_v5, :ea_v30, :ea_v32, :ea_v66]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,nothing,regressionDataset[:, :clust_id])
results3 = inference(ols) # This is assuming homoskedacity...
results3.β
results3.t

# MODEL 4 [Check]
#---------------------------------------------------------
x_name = [:FC, :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm, :ea_v5, :ea_v30, :ea_v32, :ea_v66, :logpop, :logpopsq, :popadultshare, :casino]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,nothing,regressionDataset[:, :clust_id])
results4 = inference(ols) # This is assuming homoskedacity...
results4.β
results4.t

# MODEL 5 [Done]
#---------------------------------------------------------
#Get array of fixed effects:
FixedEffects = reduce(hcat, [regressionDataset[:,:statenumber].==fe for fe in unique(regressionDataset[:,:statenumber])])
x_name = [:FC, :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm, :ea_v5, :ea_v30, :ea_v32, :ea_v66, :logpop, :logpopsq, :popadultshare, :casino]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,FixedEffects,regressionDataset[:, :clust_id])
results5 = inference(ols) # This is assuming homoskedacity...
results5.β
results5.t

# Compile results from different ols models:
table3 = [
results1.β[1] results2.β[1] results3.β[1] results4.β[1] results5.β[1]
results1.t[1] results2.t[1] results3.t[1] results4.t[1] results5.t[1]
results1.β[2] results2.β[2] results3.β[2] results4.β[2] results5.β[2]
results1.t[2] results2.t[2] results3.t[2] results4.t[2] results5.t[2]
]


# function tsls_regression(y, d, z, x)
 
#     if cov = true    
#         z = hcat(1,z,x)
#         d = hcat(1,d,x)
#     end
    
#     zᵀz = transpose(z)*z
#     zᵀd = transpose(z)*d
#     Π = inv(zᵀz) * zᵀd 
#     pred = transpose(Π)*transpose(z)
#     β = inv(pred*d)*pred*y

#     return beta

# end





