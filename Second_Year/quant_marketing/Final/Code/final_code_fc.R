

set_parameters = function(rho=.9, 
                          alpha=.2, 
                          gamma=.5, 
                          beta=.998, 
                          tau=1, 
                          sigma_nu=1, 
                          sigma_tau=1, 
                          sigma_0=1, 
                          mu_0=5,
                          vartheta = array(c(.5, .5), c(1,2))
                          ){
  
  #Set parameters for the model:

  parameters = list(rho = rho,
                    alpha = alpha, 
                    gamma = gamma, 
                    beta = beta, 
                    tau = tau, 
                    sigma_nu = sigma_nu,
                    sigma_tau = sigma_tau,
                    sigma_0 = sigma_0,
                    mu_0 = mu_0,
                    vartheta = vartheta
                    )
                    
  return(parameters)
}


params = set_parameters()

generate_fake_data = function(params){
  
  nObs = 100

  nProducts = 2
  
  P = array(runif(nProducts), c(1,nProducts))
  
  nu_jt = array(rnorm(nObs), c(nObs, 2))
  
  xi_jt = params$vartheta + nu_jt
    
}


update_quality_mean = function(){
  
  
  
}

update_quality_variance = function(){
  
  
}


choice_utility = function(){
  
  
}


softmax = function(){
  
  
}

value_function_iteration = function(){
  
  
}

