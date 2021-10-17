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
    #   xᵀx = inv(x' * x)
    #   covar = sum(u.^2) * xᵀx
    #   covar = covar .* (1 / (N - K)) # dof adjustment

    # Get standard errors, t-statistics, and p-values
    # se = sqrt.(covar[diagind(covar)])
    se = se_cluster(fit)
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

forcedcoexistence_webfinal = DataFrame(load("11423_Data_and_Programs/forcedcoexistence_webfinal.dta"));

regressionDataset = subset(forcedcoexistence_webfinal, :year => ByRow(year -> year ==(2000)))

regressionDataset[:,:clust_id] = string.(lpad.(regressionDataset[:,:eaid],3,'0')) .* string.(lpad.(regressionDataset[:,:statenumber],2,'0'))

regressionDataset[:,:instrument_precious] = regressionDataset[:,:instrument_gold] .+ regressionDataset[:,:intrument_silver]

regressionDataset[:,:wprec_enviro] = regressionDataset[:,:wgold_enviro]+regressionDataset[:,:wsilver_enviro] 


println(sum(regressionDataset[:,:instrument_precious].==0.0), " number of 0 in instrument")

# Replicate Table 3:
# Define dependent variable:
y_name = [:logpcinc]

# PANEL A:
#
# MODEL 1
#---------------------------------------------------------
# Define `independent' variables
x_name = [:FC, :HC]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,nothing,regressionDataset[:, :clust_id])
# Get standard errors:
results1a = inference(ols)

# MODEL 2 [Check]
#---------------------------------------------------------
# Define `independent' variables
x_name = [:FC, :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,nothing,regressionDataset[:, :clust_id])
results2a = inference(ols) # This is assuming homoskedacity...

# MODEL 3 [Check]
#---------------------------------------------------------
# Define `independent' variables
x_name = [:FC, :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm, :ea_v5, :ea_v30, :ea_v32, :ea_v66]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,nothing,regressionDataset[:, :clust_id])
results3a = inference(ols) # This is assuming homoskedacity...

# MODEL 4 [Check]
#---------------------------------------------------------
x_name = [:FC, :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm, :ea_v5, :ea_v30, :ea_v32, :ea_v66, :logpop, :logpopsq, :popadultshare, :casino]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,nothing,regressionDataset[:, :clust_id])
results4a = inference(ols) # This is assuming homoskedacity...

# MODEL 5 [Done]
#---------------------------------------------------------
#Get array of fixed effects:
FixedEffects = reduce(hcat, [regressionDataset[:,:statenumber].==fe for fe in unique(regressionDataset[:,:statenumber])])
x_name = [:FC, :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm, :ea_v5, :ea_v30, :ea_v32, :ea_v66, :logpop, :logpopsq, :popadultshare, :casino]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,FixedEffects,regressionDataset[:, :clust_id])
results5a = inference(ols) # This is assuming homoskedacity...


# PANEL B:
#
# MODEL 1
#---------------------------------------------------------
# Define `independent' variable
x_name = [:FC]
# Define tribe fixed effects: 
FixedEffects = reduce(hcat, [regressionDataset[:,:eaid].==fe for fe in unique(regressionDataset[:,:eaid])])

# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,FixedEffects,regressionDataset[:, :clust_id])
# Get standard errors:
results1b = inference(ols) # This is assuming homoskedacity...

# MODEL 2 
#---------------------------------------------------------
# Define `independent' variables
x_name = [:FC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,FixedEffects,regressionDataset[:, :clust_id])
results2b = inference(ols) # This is assuming homoskedacity...

# MODEL 3
#---------------------------------------------------------
# Define `independent' variables
x_name = [:FC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm, :ea_v5, :ea_v30, :ea_v66] # removed ea_v32 due to colinearity
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,FixedEffects[:,2:end],regressionDataset[:, :clust_id])
results3b = inference(ols) # This is assuming homoskedacity...


# MODEL 4 
#---------------------------------------------------------
x_name = [:FC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm, :ea_v5, :ea_v30, :ea_v66, :logpop, :logpopsq, :popadultshare, :casino]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,FixedEffects[:,2:end],regressionDataset[:, :clust_id])
results4b = inference(ols) # This is assuming homoskedacity...

# MODEL 5 
#---------------------------------------------------------
#Get array of fixed effects:
FixedEffectsState = reduce(hcat, [regressionDataset[:,:statenumber].==fe for fe in unique(regressionDataset[:,:statenumber])])
FixedEffects = hcat(FixedEffects, FixedEffectsState)
x_name = [:FC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness,  :logresarea_sqkm, :ea_v5, :ea_v30, :ea_v66, :logpop, :logpopsq, :popadultshare, :casino]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name);
# Estimate:
ols = olsRegression(x,y,FixedEffects,regressionDataset[:, :clust_id])
results5b = inference(ols) # This is assuming homoskedacity...


# Compile results from different ols models:
table3 = [
results1a.β[1] results2a.β[1] results3a.β[1] results4a.β[1] results5a.β[1]
results1a.t[1] results2a.t[1] results3a.t[1] results4a.t[1] results5a.t[1]
results1a.β[2] results2a.β[2] results3a.β[2] results4a.β[2] results5a.β[2]
results1a.t[2] results2a.t[2] results3a.t[2] results4a.t[2] results5a.t[2]
results1a.r    results2a.r    results3a.r    results4a.r    results5a.r

results1b.β[1] results2b.β[1] results3b.β[1] results4b.β[1] results5b.β[1]
results1b.t[1] results2b.t[1] results3b.t[1] results4b.t[1] results5b.t[1]
results1b.r    results2b.r    results3b.r    results4b.r    results5b.r
]

CSV.write("table3.csv", Tables.table(table3), writeheader=false) 




# Testing probit:
#----------------------------

y_name = [:FC]
x_name = [:HC, :instrument_gold, :intrument_silver]
d_name = [:HC]
# Select, drop missings and convert variables to arrays
y, x = select_variables(regressionDataset, y_name, x_name, d_name);


probit = probitModel(x,y)
probit.θ
probit.μ
probit.σ


d = y


# Testing propensityScoreMatching:
#---------------------------------

function psm_table(kn)

    # Panel A: Two instruments...
    # Subpanel 1
    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_gold, :intrument_silver , :HC];
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_a1 = propensityScoreMatching(x, y, d, kn)


    # Subpanel 2
    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_gold, :intrument_silver , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm];
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_a2 = propensityScoreMatching(x, y, d, kn)

    # Subpanel 3
    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_gold, :intrument_silver , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66];
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_a3 = propensityScoreMatching(x, y, d, kn)

    # Subpanel 4
    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_gold, :intrument_silver , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino];
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_a4 = propensityScoreMatching(x, y, d, kn)

    # Subpanel 5
    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_gold, :intrument_silver , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino];
    FixedEffects = reduce(hcat, [regressionDataset[:,:statenumber].==fe for fe in unique(regressionDataset[:,:statenumber])]);
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_a5 = propensityScoreMatching(hcat(x, FixedEffects[:,1:end-1]), y, d, kn)


    # Subpanel 6
    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_gold, :intrument_silver , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino, :removal, :homelandruggedness];
    FixedEffects = reduce(hcat, [regressionDataset[:,:statenumber].==fe for fe in unique(regressionDataset[:,:statenumber])]);
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_a6 = propensityScoreMatching(hcat(x, FixedEffects[:,1:end-1]), y, d, kn)



    # Panel B: One instrument...

    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_precious , :HC];
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_b1 = propensityScoreMatching(x, y, d, kn)


    # Subpanel 2
    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_precious , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm];
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_b2 = propensityScoreMatching(x, y, d, kn)

    # Subpanel 3
    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_precious , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66];
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_b3 = propensityScoreMatching(x, y, d, kn)

    # Subpanel 4
    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_precious , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino];
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_b4 = propensityScoreMatching(x, y, d, kn)

    # Subpanel 5
    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_precious , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino];
    FixedEffects = reduce(hcat, [regressionDataset[:,:statenumber].==fe for fe in unique(regressionDataset[:,:statenumber])]);
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_b5 = propensityScoreMatching(hcat(x, FixedEffects[:,1:end-1]), y, d, kn)


    # Subpanel 6
    y_name = [:logpcinc];
    d_name = [:FC];
    x_name = [:instrument_precious , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino, :removal, :homelandruggedness];
    FixedEffects = reduce(hcat, [regressionDataset[:,:statenumber].==fe for fe in unique(regressionDataset[:,:statenumber])]);
    y, x, d = select_variables(regressionDataset, y_name, x_name, d_name);
    psm_b6 = propensityScoreMatching(hcat(x, FixedEffects[:,1:end-1]), y, d, kn)


    # Compile results from different psm models:
    table_psm = [
        psm_a1[1] psm_a2[1] psm_a3[1] psm_a4[1] psm_a5[1] psm_a6[1]
        psm_a1[2] psm_a2[2] psm_a3[2] psm_a4[2] psm_a5[2] psm_a6[2]

        psm_b1[1] psm_b2[1] psm_b3[1] psm_b4[1] psm_b5[1] psm_b6[1]
        psm_b1[2] psm_b2[2] psm_b3[2] psm_b4[2] psm_b5[2] psm_b6[2]
    ]

    return table_psm

end

table_psm_knn3= psm_table(3)
table_psm_knn4= psm_table(4)
table_psm_knn5= psm_table(5)
table_psm_knn6= psm_table(6)


CSV.write("table_psm_knn3.csv", Tables.table(round.(table_psm_knn3,digits=3)), writeheader=false) 
CSV.write("table_psm_knn4.csv", Tables.table(round.(table_psm_knn4,digits=3)), writeheader=false) 
CSV.write("table_psm_knn5.csv", Tables.table(round.(table_psm_knn5,digits=3)), writeheader=false) 
CSV.write("table_psm_knn6.csv", Tables.table(round.(table_psm_knn6,digits=3)), writeheader=false) 





# TSLS:
#---------------------------------------------------------

# PANEL A
y_name = [:logpcinc]; d_name = [:FC];
x_name = [:instrument_gold, :intrument_silver , :HC];
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1:2]; x = H[:,3:end];
β_iv_a1 = tsls_regression(y,d,z,x) 

x_name = [:instrument_gold, :intrument_silver , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm];
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1:2]; x = H[:,3:end];
β_iv_a2 = tsls_regression(y,d,z,x) 

x_name = [:instrument_gold, :intrument_silver , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66];
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1:2]; x = H[:,3:end];
β_iv_a3 = tsls_regression(y,d,z,x) 

x_name = [:instrument_gold, :intrument_silver , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino];
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1:2]; x = H[:,3:end];
β_iv_a4 = tsls_regression(y,d,z,x) 

x_name = [:instrument_gold, :intrument_silver , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino];
FixedEffects = reduce(hcat, [regressionDataset[:,:statenumber].==fe for fe in unique(regressionDataset[:,:statenumber])]);
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1:2]; x = H[:,3:end];
β_iv_a5 = tsls_regression(y,d,z,hcat(x, FixedEffects[:,1:end-1])) 

x_name = [:instrument_gold, :intrument_silver , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino, :removal, :homelandruggedness, :wprec_enviro];
FixedEffects = reduce(hcat, [regressionDataset[:,:statenumber].==fe for fe in unique(regressionDataset[:,:statenumber])]);
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1:2]; x = H[:,3:end];
β_iv_a6 = tsls_regression(y,d,z,hcat(x, FixedEffects[:,1:end-1])) 



# PANEL B
y_name = [:logpcinc]; d_name = [:FC];
x_name = [:instrument_precious , :HC];
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1]; x = H[:,2:end];
β_iv_b1 = tsls_regression(y,d,z,x) 

x_name = [:instrument_precious , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm];
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1]; x = H[:,2:end];
β_iv_b2 = tsls_regression(y,d,z,x) 

x_name = [:instrument_precious , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66];
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1]; x = H[:,2:end];
β_iv_b3 = tsls_regression(y,d,z,x) 

x_name = [:instrument_precious , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino];
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1]; x = H[:,2:end];
β_iv_b4 = tsls_regression(y,d,z,x) 

x_name = [:instrument_precious , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino];
FixedEffects = reduce(hcat, [regressionDataset[:,:statenumber].==fe for fe in unique(regressionDataset[:,:statenumber])]);
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1]; x = H[:,2:end];
β_iv_b5 = tsls_regression(y,d,z,hcat(x, FixedEffects[:,1:end-1])) 

x_name = [:instrument_precious , :HC, :logpcinc_co, :logunempl_co, :logdist, :logruggedness, :logresarea_sqkm, :ea_v5,  :ea_v30, :ea_v32, :ea_v66, :logpop,  :popadultshare, :casino, :removal, :homelandruggedness, :wprec_enviro];
FixedEffects = reduce(hcat, [regressionDataset[:,:statenumber].==fe for fe in unique(regressionDataset[:,:statenumber])]);
y, H, d = select_variables(regressionDataset, y_name, x_name, d_name);
z = H[:,1]; x = H[:,2:end];
β_iv_b6 = tsls_regression(y,d,z,hcat(x, FixedEffects[:,1:end-1])) 


table5 = [
    β_iv_a1[1] β_iv_a2[1] β_iv_a3[1] β_iv_a4[1] β_iv_a5[1] β_iv_a6[1]
    β_iv_b1[1] β_iv_b2[1] β_iv_b3[1] β_iv_b4[1] β_iv_b5[1] β_iv_b6[1]
]

CSV.write("table5.csv", Tables.table(round.(table5,digits=3)), writeheader=false) 

