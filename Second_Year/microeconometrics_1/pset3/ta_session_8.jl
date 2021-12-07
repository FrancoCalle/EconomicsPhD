using MyMethods
using LinearAlgebra, Statistics
using JuMP, GLPK # for linear optimization
using Plots, LaTeXStrings # for plots



# Simulate from DGP

# Set number of draws from U for integral computation
nMC = 10000

# Define Z and p as in the paper
Z = collect(1:4)
p = [0.12, 0.29, 0.48, 0.78]

# Draw from U for integral computation
u = collect(1:nMC)./(nMC+1) #rand(N)
ind_1 = reduce(hcat, map(u -> u .<= p,u))' # 1{u <= P(z)}

# MTRs
m0 = 0.9 .- 1.1u .+ 0.3u.^2;
m1 = 0.35 .- 0.3u .- 0.05u.^2;



# Recreate Figure 1

# Initialize plot
pyplot(size=(500,300), leg=true);

# Plot parameters
plt=plot(u, m1, line = (:line, :line, 1, 1, :orange), label=L"$m_1(u)$")
plot!(u, m0, line = (:line, :line, 1, 1, :green), label=L"$m_0(u)$")
plot!(u, m1 - m0, line = (:line, :line, 1, 1, :blue), label=L"$m_1(u) - m_0(u)$")

# Additional formatting
plot!(ylim=[-1,1], yticks = -1:0.2:1, xticks = 0:0.1:1, 
    legend=(0.69,0.675), framestyle = :box, grid=false, 
    background_color_legend = nothing)
ylabel!("Parameter value")
xlabel!("Unobserved heterogeneity in treatment choice "*L"$(u)$")



# Get ATT weights
w1_ATT = ind_1 ./ mean(p)
w0_ATT = -w1_ATT;

# Calculate the implied ATT
b_ATT = mean(w1_ATT .* m1) + mean(w0_ATT .* m0)
print(round(b_ATT, digits = 3))



# Compute β_s

# IV slope
Cov_DZ = 0.25 *(sum(Z.*p)) - mean(p)*mean(Z) # Define population moments
s_IV = ((Z .- mean(Z))./Cov_DZ)' # s(d,z)
w1_IV = s_IV .* ind_1 # weight for d = 1
w0_IV = s_IV .* (1 .- ind_1) # weight for d = 0
b_IV = mean(w1_IV .* m1) + mean(w0_IV .* m0); # implied IV slope in DGP

# TSLS slope
ZD = hcat(ones(4)*0.25, 0.25.*p)' # Define population moments
inv_ZZ = I.*4
FS = inv_ZZ * ZD # first stage
s_TSLS = (inv(FS * ZD') * ZD * inv_ZZ)[2,:]' # s(d,z)
w1_TSLS = s_TSLS .* ind_1 # weight for d = 1
w0_TSLS = s_TSLS .* (1 .- ind_1) # weight for d = 0
b_TSLS = mean(w1_TSLS .* m1) + mean(w0_TSLS .* m0); # implied IV slope in DGP

# Print β_s
print((b_IV=round(b_IV, digits = 3), b_TSLS=round(b_TSLS, digits = 3)))



# function to get gammas
function get_gamma(B, w1, w0)
    Bw1 = B .* mapslices(mean, w1, dims=2)
    Bw0 = B .* mapslices(mean, w0, dims=2)
    return mapslices(mean, Bw1, dims=1), mapslices(mean, Bw0, dims=1)
end



# Function to obtain lower/upper bound for ATT given IV and TSLS slope
function get_bound(b_IV, b_TSLS,
                    gamma_1IV, gamma_0IV,
                    gamma_1TSLS, gamma_0TSLS,
                    gamma_1ATT, gamma_0ATT;
                    sense="Max", decreasing=false)
                    
    # Data parameters
    K = length(gamma_1IV) -1

    # initialize model
    m = Model(GLPK.Optimizer)

    # initialize variables
    @variable(m, 1 >= theta[1:(K+1), 1:2] >= 0) # bounded by 0 and 1

    # set constraints
    @constraint(m, IV, sum(theta[:,1].*gamma_1IV') +
        sum(theta[:,2].*gamma_0IV')  == b_IV)
    @constraint(m, TSLS, sum(theta[:,1].*gamma_1TSLS') +
        sum(theta[:,2].*gamma_0TSLS') == b_TSLS)

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


# Compute bounds with flexible parametric specification

# Compute bounds for different polynomial degrees
p_bounds = zeros(2, 2, 19)
for k in 1:19
    # Get Bernstein poly basis terms and compute the gamma
    Bernstein = MyMethods.get_basis(u, "Bernstein", k, nothing)
    gamma_1IV, gamma_0IV = get_gamma(Bernstein, w1_IV, w0_IV)
    gamma_1TSLS, gamma_0TSLS = get_gamma(Bernstein, w1_TSLS, w0_TSLS)
    gamma_1ATT, gamma_0ATT = get_gamma(Bernstein, w1_ATT, w0_ATT)

    # Compute the parametric bounds
    for dec in (false,true)
        # lower bound
        p_bounds[dec+1, 1, k] = get_bound(b_IV, b_TSLS,
                gamma_1IV, gamma_0IV,
                gamma_1TSLS, gamma_0TSLS,
                gamma_1ATT, gamma_0ATT,
                sense = "Min", decreasing = dec)

        # upper bound
        p_bounds[dec+1, 2, k] = get_bound(b_IV, b_TSLS,
                gamma_1IV, gamma_0IV,
                gamma_1TSLS, gamma_0TSLS,
                gamma_1ATT, gamma_0ATT,
                sense = "Max", decreasing = dec)
    end
end



# Get Bernstein poly basis terms and compute the gamma
k = 4
Bernstein = MyMethods.get_basis(u, "Bernstein", k, nothing)
gamma_1IV, gamma_0IV = get_gamma(Bernstein, w1_IV[:,1], w0_IV[:,1])
gamma_1TSLS, gamma_0TSLS = get_gamma(Bernstein, w1_TSLS, w0_TSLS)
gamma_1ATT, gamma_0ATT = get_gamma(Bernstein, w1_ATT, w0_ATT)

# Compute the parametric bounds
p_bounds = zeros(2, 2, 19)

for dec in (false,true)
    # lower bound
    p_bounds[dec+1, 1, k] = get_bound(b_IV, b_TSLS,
            gamma_1IV, gamma_0IV,
            gamma_1TSLS, gamma_0TSLS,
            gamma_1ATT, gamma_0ATT,
            sense = "Min", decreasing = dec)

    # upper bound
    p_bounds[dec+1, 2, k] = get_bound(b_IV, b_TSLS,
            gamma_1IV, gamma_0IV,
            gamma_1TSLS, gamma_0TSLS,
            gamma_1ATT, gamma_0ATT,
            sense = "Max", decreasing = dec)
end

p_bounds[:,:,5]



# Compute nonparametric bounds

# Get constant splines (indicator specificaton) and compute the gamma
CSplines = MyMethods.get_basis(u, "CSplines2", 3, p)[:,1:4] 
gamma_1IV, gamma_0IV = get_gamma(CSplines, w1_IV, w0_IV)
gamma_1TSLS, gamma_0TSLS = get_gamma(CSplines, w1_TSLS, w0_TSLS)
gamma_1ATT, gamma_0ATT = get_gamma(CSplines, w1_ATT, w0_ATT)

# Compute the nonparametric bounds
np_bounds = zeros(2,2)
for dec in (false,true)
    # lower bound
    np_bounds[dec+1, 1] = get_bound(b_IV, b_TSLS,
            gamma_1IV, gamma_0IV,
            gamma_1TSLS, gamma_0TSLS,
            gamma_1ATT, gamma_0ATT,
            sense = "Min", decreasing = dec)
    
    # upper bound
    np_bounds[dec+1, 2] = get_bound(b_IV, b_TSLS,
            gamma_1IV, gamma_0IV,
            gamma_1TSLS, gamma_0TSLS,
            gamma_1ATT, gamma_0ATT,
            sense = "Max", decreasing = dec)
end

# cast np_bounds to dimension of p_bounds
np_bounds = reshape(repeat(np_bounds, 19)', (2, 2, 19));




# Plot the bounds

# Initialize plot
pyplot(size=(500,300), leg=true);
_x = collect(1:19)

# ATT 
plt=plot(_x,ones(19)*b_ATT, line = (:line, :dot, 0.8, 1, :black), label="")
plot!(_x,ones(19)*b_ATT, seriestype = :scatter, markershape=:star5, markersize=4, color=:black, label="")
annotate!((19.24,-0.325, text("ATT", :black, :right, 8)))

# parametric bounds
plot!(_x, p_bounds[1,1,:], line = (:line, :line, 0.8, 1, :blue), label = "Polynomial")
plot!(_x, p_bounds[1,1,:], seriestype = :scatter, markershape=:circle, markersize=4, color=:blue, label="")
plot!(_x, p_bounds[1,2,:], line = (:line, :line, 0.8, 1, :blue), label="")
plot!(_x, p_bounds[1,2,:], seriestype = :scatter, markershape=:circle, markersize=4, color=:blue, label="")

# nonparametric bounds
plot!(_x, np_bounds[1,1,:], line = (:line, :dot, 0.8, 1, :blue), label="Nonparametric")
plot!(_x, np_bounds[1,1,:], seriestype = :scatter, markershape=:circle, markersize=4, color=:blue, label="")
plot!(_x, np_bounds[2,1,:], line = (:line, :dot, 0.8, 1, :blue), label="")
plot!(_x, np_bounds[2,1,:], seriestype = :scatter, markershape=:circle, markersize=4, color=:blue, label="")

# parametric bounds w/ decreasing MTRs
plot!(_x, p_bounds[2,1,:], line = (:line, :line, 0.8, 1, :orange), label="Polynomial and decreasing")
plot!(_x, p_bounds[2,1,:], seriestype = :scatter, markershape=:rect, markersize=4, color=:orange, label="")
plot!(_x, p_bounds[2,2,:], line = (:line, :line, 0.8, 1, :orange), label="")
plot!(_x, p_bounds[2,2,:], seriestype = :scatter, markershape=:rect, markersize=4, color=:orange, label="")

# nonparametric bounds w/ decreasing MTRs
plot!(_x, np_bounds[1,2,:], line = (:line, :dot, 0.8, 1, :orange), label="Nonparametric and decreasing")
plot!(_x, np_bounds[1,2,:], seriestype = :scatter, markershape=:rect, markersize=4, color=:orange, label="")
plot!(_x, np_bounds[2,2,:], line = (:line, :dot, 0.8, 1, :orange), label="")
plot!(_x, np_bounds[2,2,:], seriestype = :scatter, markershape=:rect, markersize=4, color=:orange, label="")

# Additional formatting
plot!(ylim=[-1,0.4], yticks = -1:0.2:0.4, xticks = 1:1:19, 
    legend=(0,0), framestyle = :box, grid=false, 
    background_color_legend = nothing)
ylabel!("Upper and lower bounds")
xlabel!("Polynomial degree "*L"$(K_0=K_1\equiv K)$")


