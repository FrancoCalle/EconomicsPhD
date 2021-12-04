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



# Compute \beta_s

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







