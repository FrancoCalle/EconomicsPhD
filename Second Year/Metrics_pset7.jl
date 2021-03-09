using Random, Distributions, LinearAlgebra, Optim, NLSolversBase
Random.seed!(11121996)

function generate_errors(N, K)
    normal_noise = MvNormal(zeros(K), ones(K))
    errors = rand(normal_noise, N)'
    return errors
end

function generate_modelled_data(N::Int,
                                K::Int,
                                A::Array,
                                F::Array,
                                x_0::Int,
                                y_0::Int)

    errors = generate_errors(N, K)
    y = Array{Float64, 2}(undef, N+1, 2)
    x = Array{Float64, 2}(undef, N+1, 2)

    y[N, 1] = y_0
    y[N, 2] = x_0
    for n = 1:(N-1)
        y[N-n, 1:2] = A*y[(N+1)-n, 1:2] + F * errors[(N+1)-n, 1:2]
        x[N-n, 1:2] = y[(N+1)-n, 1:2]
    end
    # Cutting off the lagged bits we don't want
    return y[2:(N-1), :], x[2:(N-1), :]
end


function multivariate_normal_distribution(z, A, Ω, K)
    Ψ = (2*pi)^(-K/2)*det(Ω)^(-1/2)*exp.(-0.5*sum(z*inv(Ω).*z,dims=2))
    return Ψ
end

# sum(z*inv(Ω).*z,dims=2)

function log_likelihood(data, A, F)

    # splitting y and x out
    y = data[1]
    x = data[2]
    # Creating dimensions
    N = size(y)[1]
    K = size(y)[2]
    # Creating variance matrix
    Ω = F * F'
    #Placeholder for Log_Li variable:
    z = y' .- A*x'

    Li = multivariate_normal_distribution(z', A, Ω, K)

    log_Li = log.(Li)

    return log_Li
end

function objective_function(data, A, F)
    log_Li = log_likelihood(data, A, F)
    return sum(log_Li)/size(data[1])[1]
end

function alpha_anon_LL(α_1)
    A_hat = [α_1 -0.2; 0 0.8]
    anon_LL = objective_function(sim_data, A_hat, generative_F_matrix)
    return -anon_LL
end

function all_anon_LL(param_array)
    A_hat = [param_array[1] param_array[2];
            0 param_array[3]]
    anon_LL = objective_function(sim_data, A_hat, generative_F_matrix)
    return -anon_LL
end



#Find Score Function for each yt
function find_Score_i(y,x,params)

    K = 2
    F = generative_F_matrix
    Ω = F * F'

    function all_anon_Log_Li(param_array)


        #Placeholder for Log_Li variable:
        A_hat = [param_array[1] param_array[2];
                0 param_array[3]]

        z = y .- (A_hat*x)


        Ψ = (2*pi)^(-K/2)*det(Ω)^(-1/2)*exp(-0.5*z'*inv(Ω)*z)

        return log(Ψ)
    end

    od = OnceDifferentiable(vars -> all_anon_Log_Li(vars[1:3]), params; autodiff=:forward)
    grad = zeros(3)
    return od.df(grad, params)

end


## Generate Fake Data

N=1000
generative_A_matrix = [0.9 -0.2; 0 0.8]
generative_F_matrix = [1 0.5; 0 1]

sim_data = generate_modelled_data(N,
                                  2,
                                  generative_A_matrix,
                                  generative_F_matrix,
                                  0,
                                  0)


## trivial

check_a = log_likelihood(sim_data, generative_A_matrix, generative_F_matrix)
check_b = objective_function(sim_data, generative_A_matrix, generative_F_matrix)

check_a = objective_function(sim_data, ones(2,2), generative_F_matrix)
check_b = objective_function(sim_data, generative_A_matrix, generative_F_matrix)

if (check_b > check_a)
    println("Success, so far")
end

##

init_params = 0.0
init_params_all = [0.0, 0.0, 0.0]

opt_alpha_1 = optimize(alpha_anon_LL, [init_params], LBFGS(),  Optim.Options(store_trace=true, extended_trace=true))
opt_all = optimize(all_anon_LL, init_params_all, LBFGS())

alpha_1_params = Optim.minimizer(opt_alpha_1)

all_params = Optim.minimizer(opt_all)

## Compute score and gradients and information under the nuisance parameters:

N = size(sim_data[1])[1]

d_LogLi = Array{Float64, 2}(undef, N, 3)
for ii =1:N
    d_LogLi[ii,:] = find_Score_i(sim_data[1][ii,:] , sim_data[2][ii,:],all_params)
end

s_prime = d_LogLi

ΔS = s_prime

ΔSα = ΔS[:,1]
ΔSϑ = ΔS[:,2:3]

β = ΔSϑ\ΔSα
U_hat = ΔSα .- ΔSϑ*β


sum(sim_data[1][:,1].^2)/N

Information_alpha_1 = mean(ΔSα.^2)
Information_nuisance = mean(U_hat.^2)
