using Pkg
Pkg.activate(".") # Create new environment in this folder

# first time you need to install dependencies
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add(["DataFrames","DataFramesMeta","Chain"])
Pkg.add("Plots")
Pkg.add("CategoricalArrays")

# past the first time, you only need to instanciate the current folder
Pkg.instantiate() # Updates packages given .toml file

#----------------------------------------------------

using Distributions
using LinearAlgebra
using StatsBase
using DataFrames
using Plots
using CategoricalArrays

#----------------------------------------------------

α_sd = 1
ψ_sd = 1


csort = 0.5 # Sorting effect
csig  = 0.5 # Cross-sectional standard deviation
w_sigma = 0.2


function data_generating_process(α_sd,
                                ψ_sd,
                                csort, 
                                csig,
                                w_sigma)

    cnetw = 0.2 # Network effect

    nk = 30
    nl = 10

    # Let's assume moving probability is fixed
    λ = 0.1

    # approximate each distribution with some points of support
    ψ = quantile.(Normal(), (1:nk) / (nk + 1)) * α_sd
    α = quantile.(Normal(), (1:nl) / (nl + 1)) * ψ_sd

    # Let's create type-specific transition matrices
    # We are going to use joint normals centered on different values
    G = zeros(nl, nk, nk)
    for l in 1:nl, k in 1:nk
        G[l, k, :] = pdf( Normal(0, csig), ψ .- cnetw * ψ[k] .- csort * α[l])
        G[l, k, :] = G[l, k, :] ./ sum(G[l, k, :])
    end

    # We then solve for the stationary distribution over psis for each alpha value
    # We apply a crude fixed point approach
    H = ones(nl, nk) ./ nk
    for l in 1:nl
        M = transpose(G[l, :, :])
        for i in 1:100
            H[l, :] = M * H[l, :]
        end
    end

    # p1 = plot(G[1, :, :], xlabel="Previous Firm", ylabel="Next Firm", zlabel="G[1, :, :]", st=:wireframe)
    # p2 = plot(G[nl, :, :], xlabel="Previous Firm", ylabel="Next Firm", zlabel="G[nl, :, :]", st=:wireframe, right_margin = 10Plots.mm) # right_margin makes sure the figure isn't cut off on the right
    # plot(p1, p2, layout = (1, 2), size=[600,300])

    #--------------------------------------------------------------

    nt = 10
    ni = 10000

    # We simulate a balanced panel
    ll = zeros(Int64, ni, nt) # Worker type
    kk = zeros(Int64, ni, nt) # Firm type
    spellcount = zeros(Int64, ni, nt) # Employment spell

    for i in 1:ni
        
        # We draw the worker type
        l = rand(1:nl)
        ll[i,:] .= l
        
        # At time 1, we draw from H
        kk[i,1] = sample(1:nk, Weights(H[l, :]))
        
        for t in 2:nt
            if rand() < λ
                kk[i,t] = sample(1:nk, Weights(G[l, kk[i,t-1], :]))
                spellcount[i,t] = spellcount[i,t-1] + 1
            else
                kk[i,t] = kk[i,t-1]
                spellcount[i,t] = spellcount[i,t-1]
            end
        end
        
    end

    #------------------------------------------------------------

    firms_per_type = 15
    jj = zeros(Int64, ni, nt) # Firm identifiers

    draw_firm_from_type(k) = sample(1:firms_per_type) + (k - 1) * firms_per_type

    for i in 1:ni
        
        # extract firm type
        k = kk[i,1]
        
        # We draw the firm (one of firms_per_type in given group)
        jj[i,1] = draw_firm_from_type(k)
        
        for t in 2:nt
            if spellcount[i,t] == spellcount[i,t-1]
                # We keep the firm the same
                jj[i,t] = jj[i,t-1]
            else
                # We draw a new firm
                k = kk[i,t]
                
                new_j = draw_firm_from_type(k)            
                # Make sure the new firm is actually new
                while new_j == jj[i,t-1]
                    new_j = draw_firm_from_type(k)
                end
                
                jj[i,t] = new_j
            end
        end
    end
    # Make sure firm ids are contiguous
    contiguous_ids = Dict( unique(jj) .=> 1:length(unique(jj))  )
    jj .= getindex.(Ref(contiguous_ids),jj);

    #----------------------------------------------------------------------------------------

    ii = repeat(1:ni,1,nt)
    tt = repeat((1:nt)',ni,1)
    df = DataFrame(i=ii[:], j=jj[:], l=ll[:], k=kk[:], α=α[ll[:]], ψ=ψ[kk[:]], t=tt[:], spell=spellcount[:]);

    return df, G, H

end



using Chain
using DataFramesMeta

#Number of workers in the firm:

crossSectionMean = @chain df begin
    groupby([:j, :t])
    combine(nrow => :count)
    groupby(:j)
    combine(:count => mean)
    @aside println("On average each firm has $(round(mean(_.count_mean))) observations, averaging across time periods.")
    # first(_, 5)
end


#Number of movers at the firm, defined as people that moved into the firm:

println("Maximum number of moves by a single worker: ", maximum(df.spell))

moversDataFrame = @chain df begin
    groupby([:j, :i, :spell])
    combine(:t => mean)
    @transform!(:spell_true = :spell .> 0)
    groupby(:j)
    combine(:spell_true => sum)
    @aside println("On average each firm has $(round(mean(_.spell_true_sum))) movers over the 10 periods.")
end


#---------------------------------------------------------------------------------------------------


df[!, :lw] = df.α + df.ψ + w_sigma * rand(Normal(), size(df)[1]);
df[!, :dummy] .= 1;


#1. Mean wage within firm: 
eventStudyPanel =  @chain df begin
    groupby([:j, :t]) # Group by firm and time 
    transform(:lw => mean) # Obtain wage at that group level
    @transform!(:wage_percentile = cut(:lw_mean,4)) # Get wage_percentile at same level of agregation
    groupby([:i, :spell])
    transform(nrow => :years_at_firm) # Get number of years individual spend at a single firm
    sort([:i,:t], rev = true) # Sort over time by individual (this helps to get first 2 periods of last firm)
    groupby([:i, :spell])  # by group and spell 
    transform(:dummy .=>  cumsum => :years_at_firm_sofar) # Get the number of years spent at a firm so far
    sort([:i,:t], rev = false) # Sort over time by individual (this helps to get first 2 periods of last firm)
    groupby([:i, :spell])  # by group and spell
    transform(:dummy .=>  cumsum => :years_at_firm_sofar_inverse) # Get the number of years spent at a firm so far, but backwards.
    @aside initialFirmDataFrame = @chain _ begin
        subset(:spell => ByRow(==(1)), :years_at_firm_sofar_inverse => ByRow(<=(2))) # Generate dataframe for initial firm, keep last two years (call it initialFirmDataFrame)
    end
    subset(:spell => ByRow(==(0)), :years_at_firm_sofar => ByRow(<=(2))) # Generate dataframe for subsequent firm, keep first two years
    append!(initialFirmDataFrame) # Append both dataframes
    groupby(:i) # group by person 
    transform(nrow => :nreps) # and get number of repetitions
    subset(:nreps => ByRow(==(4))) # and get number of repetitions
    sort([:i,:t]) # sort by i and t to check if it worked out
    # Generating event time variable:
    sort([:i,:t], rev = false) # Sort over time by individual (this helps to get first 2 periods of last firm)
    groupby([:i])  # by worker and time
    transform(:dummy .=>  cumsum => :event_time) # Get the number of years spent at a firm so far, but backwards.
end

# Generate dictionary to index cut label:
percentile_cut = Dict(1:4 .=> unique(sort(eventStudyPanel.wage_percentile)))

# Get Card Heining and Kling event Figure:
# There may be a smart way to do this, here I'm brute forcing it ...
initial = 1
final = 4
function generateEventStudy(eventStudyPanel, initial, final)
    return  eventStudy =  @chain eventStudyPanel begin
                            @aside finalJob = @chain _ begin
                                subset(:wage_percentile => ByRow(==(percentile_cut[final])), :spell => ByRow(==(1)))
                            end
                            subset(:wage_percentile => ByRow(==(percentile_cut[initial])), :spell => ByRow(==(0)))
                            append!(finalJob) # Append both dataframes
                            groupby(:i)
                            transform(nrow => :nreps)
                            subset(:nreps => ByRow(==(4)))
                            groupby([:wage_percentile, :event_time])
                            combine(:lw => mean)
                            sort(:event_time)
    end
end


eventStudySwitchers = [generateEventStudy(eventStudyPanel, 1, 1).lw_mean,
                        generateEventStudy(eventStudyPanel, 1, 2).lw_mean,
                        generateEventStudy(eventStudyPanel, 1, 3).lw_mean,
                        generateEventStudy(eventStudyPanel, 1, 4).lw_mean]
eventStudySwitchers2 = [generateEventStudy(eventStudyPanel, 4, 1).lw_mean,
                        generateEventStudy(eventStudyPanel, 4, 2).lw_mean,
                        generateEventStudy(eventStudyPanel, 4, 3).lw_mean,
                        generateEventStudy(eventStudyPanel, 4, 4).lw_mean];

plot(1:4,eventStudySwitchers, label=["1 to 1" "1 to 2" "1 to 3" "1 to 4"],markershape = :square)
plot!(1:4,eventStudySwitchers2,label =  ["4 to 1 " "4 to 2" "4 to 3" "4 to 4"],markershape = :circle)
plot!(legend=:outertopright)


# Calibrating the parameters:
#----------------------------------------------------------------------------------------

# Generate function with variance decomposition:

function variance_decomposition(df)
    varianceDecomposition =  @chain df begin
        groupby([:k]) # Group by firm type 
        combine(:α => mean, :ψ => mean)
    end

    std_α = std(varianceDecomposition.α_mean)
    std_ψ = std(varianceDecomposition.ψ_mean)
    std_αψ = cov(varianceDecomposition.α_mean,varianceDecomposition.ψ_mean)

    return [std_α std_ψ std_αψ]
end

α_sd = 1
ψ_sd = 1
csort = 0.5 # Sorting effect
csig  = 0.5 # Cross-sectional standard deviation
w_sigma = 0.2

gap = 0.20
α_sd_list = 0.01:gap:1
ψ_sd_list = 0.01:gap:1
csort_list = 0.01:gap:1 # Sorting effect
csig_list  = 0.01:gap:1 # Cross-sectional standard deviation
w_sigma_list = 0.01:gap:1

flat_gridpoints(grids) = vec(collect(Iterators.product(grids...)))
grid_points = flat_gridpoints((α_sd_list, ψ_sd_list, csort_list, csig_list,w_sigma_list));

jj = 1

true_var_decomp = [0.186 0.101 0.110]
difference_outcome = []

#Apply grid search, this can be optimized for sure...
for jj in 1:length(grid_points)
    α_sd,ψ_sd,csort,csig, w_sigma = grid_points[jj]
    df,_,_ = data_generating_process(α_sd,ψ_sd,csort,csig, w_sigma);
    var_decomp = variance_decomposition(df)
    euclidean_dist = sum((true_var_decomp-var_decomp).^2)
    push!(difference_outcome, euclidean_dist)
end


# Get index of minimum euclidean distance estimation:
_, index_calibration = findmin(difference_outcome)
α_sd,ψ_sd,csort,csig, w_sigma = grid_points[index_calibration]
df, G, H  = data_generating_process(α_sd,ψ_sd,csort,csig, w_sigma);


p1 = plot(G[1, :, :], xlabel="Previous Firm", ylabel="Next Firm", zlabel="G[1, :, :]", st=:wireframe)
p2 = plot(G[nl, :, :], xlabel="Previous Firm", ylabel="Next Firm", zlabel="G[nl, :, :]", st=:wireframe, right_margin = 10Plots.mm) # right_margin makes sure the figure isn't cut off on the right
plot(p1, p2, layout = (1, 2), size=[600,300])

plot(H, xlabel="Worker", ylabel="Firm", zlabel="H", st=:wireframe)


# The output reduces cross sectional variance too much, check what happens when ↑ csig ...
df, G, H  = data_generating_process(α_sd,ψ_sd,csort, .1, w_sigma);
var_decomp = variance_decomposition(df)
euclidean_dist = sum((true_var_decomp-var_decomp).^2)
println("Euclidean distance is: $euclidean_dist")

p1 = plot(G[1, :, :], xlabel="Previous Firm", ylabel="Next Firm", zlabel="G[1, :, :]", st=:wireframe)
p2 = plot(G[nl, :, :], xlabel="Previous Firm", ylabel="Next Firm", zlabel="G[nl, :, :]", st=:wireframe, right_margin = 10Plots.mm) # right_margin makes sure the figure isn't cut off on the right
plot(p1, p2, layout = (1, 2), size=[600,300])

# It seems csig contributes importantly to covariance! Check this again and get feedback from prof.
# I'll leave it at 0.05 which is 5 times it optimal* value.


# Delete method in case of mistake
# m = @which data_generating_process(α_sd,ψ_sd,csort,csig, w_sigma)
# Base.delete_method(m)

#--------------------------------------------------------------------------------------

# Estimating two way fixed effects ... 
Pkg.add("LightGraphs")
# Pkg.add("TikzGraphs")


using LightGraphs
using TikzGraphs

# First create a matrix [NxT, nFirms, nFirms ] indicating  where each single individual is going...
nfirms = length(unique(df.j))

function individualDeterministicTransitionMatrix(df,ii)

    nfirms = length(unique(df.j))
    shift_i = @chain df begin
                        subset(:i => ByRow(.==(ii)))
                        sort(:t)
                        combine(first, groupby(_,:spell))
                    end
    
    nshifts = size(shift_i)[1]
    transitionMatrix_i = zeros(Int32, nfirms, nfirms);
    
    for ii in 1:nshifts
        if ii != nshifts
            current = shift_i.j[ii]
            next = shift_i.j[ii+1]
            transitionMatrix_i[current,next] = 1
        end
    end

    return transitionMatrix_i
end

# Add all shifts happening in the economy during the whole time
totalDeterministicShifts = zeros(nfirms, nfirms);
for ii in unique(df.i)
    totalDeterministicShifts = totalDeterministicShifts + individualDeterministicTransitionMatrix(df,ii)
end

adjacencyMatrix = (UpperTriangular(totalDeterministicShifts) + transpose(LowerTriangular(totalDeterministicShifts))).>1
adjacencyMatrix = adjacencyMatrix + transpose(adjacencyMatrix)

println("There are $(sum(totalDeterministicShifts)) job shifts across time ...")
println("There are $(sum(adjacencyMatrix)) edges between nodes ...")

simpleGraph = SimpleGraph(adjacencyMatrix);
connectedNetwork = connected_components(g)
connectedSet = connectedNetwork[1]

println("We have only $(length(connectedSet)) firms fully connected ...") # Double check we might be wrong...

df_connected = df[in(connectedSet).(df.j),:]
println("Due to unconnectedness we eliminated $(size(df)[1]-size(df_connected)[1]) observations, not much")


#--------------------------------------------------------------------------------------


