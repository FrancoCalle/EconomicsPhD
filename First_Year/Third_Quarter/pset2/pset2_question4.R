# Created by: Franco Calle

# Pset 2, Question 4

# Preamble
#---------

packages = c("tidyverse", "haven", "knitr", 
             "MatchIt", "sandwich", "stargazer", 
             "gridExtra","ggplot2", "optimization")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

root <- "C:/Users/franc/Dropbox/Franco First Year/Empirical Analysis III/problem_sets/pset2/question4"

# Reading in data and setting up control variables
incomeData <- 
  read_csv(file.path(root, "IncomeData.csv")) 

# Set parameters
lb = 18000
ub = 22000
h = (ub - lb)/2
K = lb + h
rho1 = 1
rho2 = .8

# Part A: Plot histogram of the data:
#------------------------------------

histogram_4a <- qplot(incomeData$Y, geom="histogram", binwidth = 500)+
  geom_vline(xintercept = K - h, linetype="dashed", color = "red", size=1)+ 
  geom_vline(xintercept = K + h, linetype="dashed", color = "red", size=1)+ 
  xlim(0, 120000) + 
  xlab("") +       
  ylab("")       

ggsave(histogram_4a, file = file.path(root, "4a_histogram.png"), width = 10, height = 6)

# Part B: Compute estiamtes for F-(18000) and F+(22000):
#-------------------------------------------------------

# incomeData$y_bin100 <- round(incomeData$Y/100,0)*100
incomeData$y_bin100 <- floor(incomeData$Y/100)*100  + 100

# Get $100 binned data and weights for each bin:
incomeDataBinned<-  incomeData %>%
                    count(y_bin100) %>%
                    mutate(frequency = prop.table(n))

ggplot(data=incomeDataBinned, aes(x=y_bin100, y=frequency)) +
  geom_bar(stat="identity") +
  geom_vline(xintercept = K - h, linetype="dashed", color = "red", size=1)+ 
  geom_vline(xintercept = K + h, linetype="dashed", color = "red", size=1)+ 
  xlim(0, 120000) + 
  ylim(0, 0.005) +
  xlab("")       


count_lb<- (
  incomeData %>% 
    filter(between(Y, lb-h, lb))
)

count_ub<- (
  incomeData %>% 
    filter(between(Y, ub, ub+h))
)



# Compute density lower bound f^-(18000)
H_lb<- (
    incomeDataBinned %>% 
    filter(between(y_bin100, lb-h, lb)) %>% 
    summarise(mean = sum(frequency))
    )[[1]]

f_lb <- H_lb/h

# Compute density upper bound f^-(22000)
H_ub<- (
    incomeDataBinned %>% 
    filter(between(y_bin100, ub, ub+h)) %>% 
    summarise(mean = sum(frequency))
    )[[1]]

f_ub <- H_ub/h

# Compute density of excluded range \int_lb^ub
H_K <- (
    incomeDataBinned %>% 
    filter(between(y_bin100, lb, ub)) %>% 
    summarise(bunchDensity = sum(frequency))
    )[[1]]


f_lb <- (nrow(count_lb)/nrow(incomeData))/h
(nrow(count_ub)/nrow(incomeData))/h

f_ub <- (nrow(count_ub)/nrow(incomeData))/h
(nrow(count_ub)/nrow(incomeData))/h


# Bunching probability:
PK = H_K - (f_lb + f_ub)*h

# PK = H_K - (H_lb + H_ub)


# Part C: Compute taxable income elasticity:
#------------------------------------------

# Density version 1:
density_1 <- function(bounds, density_bounds){

  eta_lb <- bounds[1] 
  eta_ub <- bounds[2] 
  phi_lb <- density_bounds[1]
  phi_ub <- density_bounds[2]

  integrand <- function(eta){
    phi <- phi_lb + (eta - eta_lb)/(eta_ub - eta_lb) * (phi_ub - phi_lb)
  }
  
  phi_int <- integrate(integrand, lower = eta_lb, upper = eta_ub)[[1]]
  
  return(phi_int)
  }

# Density version 2:
density_2 <- function(bounds, density_bounds, threshold){
  
  eta_lb <- bounds[1] 
  eta_ub <- bounds[2] 
  phi_lb <- density_bounds[1]
  phi_ub <- density_bounds[2]

  
  integrand <- function(eta){
    phi <- ifelse(eta <= threshold, 
                  phi_lb, 
                  phi_lb + (eta - threshold)/(eta_ub - threshold) * (phi_ub - phi_lb)
                  )
    }

  phi_int <- integrate(integrand, lower = eta_lb, upper = eta_ub)[[1]]
  
  return(phi_int)
}

# Density version 3:
density_3 <- function(bounds, density_bounds, threshold){

  eta_lb <- bounds[1] 
  eta_ub <- bounds[2] 
  phi_lb <- density_bounds[1]
  phi_ub <- density_bounds[2]


  integrand <- function(eta){
    phi <- ifelse(eta <= threshold, 
                  phi_lb + (eta - eta_lb)/(threshold - eta_lb) * (phi_ub - phi_lb), 
                  phi_ub)
  }
  
  phi_int <- integrate(integrand, lower = eta_lb, upper = eta_ub)[[1]]
  
  return(phi_int)
}


# Compute density over interval

preference_density <- function(bbeta, 
                               version,
                               f_lb, 
                               f_ub){
  
  # Compute eta lower and upper bound
  eta_lb <- lb * rho1 ^ (-bbeta)
  eta_ub <- ub * rho2 ^ (-bbeta)
  preference_bounds <- c(eta_lb, eta_ub)
  
  # Compute density evaluated in upper bound and lower bound
  phi_lb <- f_lb * rho1 ^ bbeta #preference_density_lb(bbeta, f_lb)
  phi_ub <- f_ub * rho2 ^ bbeta #preference_density_ub(bbeta, f_ub)
  density_bounds <- c(phi_lb, phi_ub)
  

  # Apply density version #:

  if(version == 1){
    phi_int <- density_1( 
            bounds=preference_bounds, 
            density_bounds=density_bounds
           ) #Apply density version 1
  } else if (version == 2){
    p <- .25
    t <- p*eta_lb + (1-p)*eta_ub
    phi_int <- density_2(
            preference_bounds, 
            density_bounds,           
            t
           ) #Apply density version 2
  } else if (version == 3){
    p <- .75
    t <- p*eta_lb + (1-p)*eta_ub
    phi_int <- density_3(
            preference_bounds, 
            density_bounds,           
            t) #Apply density version 3
  }
  
    
  return(phi_int)
}



#Optimize problem:


# Version 1:
loss_v1 = function(bbeta, vv, f_lb, f_ub, PK) {
  (PK-preference_density(bbeta, vv, f_lb,f_ub))^2
}

ov1 = optim(0.2, 
            loss_v1, 
            vv=1, #Function parameters
            f_lb=f_lb, 
            f_ub=f_ub, 
            PK=PK, #Optimization method
            method="Brent",
            lower = 0, 
            upper = 5)

#Version 2:
loss_v2 = function(bbeta, vv, f_lb, f_ub, PK) {
  (PK-preference_density(bbeta, vv, f_lb,f_ub))^2
}

ov2 = optim(0.2, loss_v2, vv=2, f_lb=f_lb, f_ub=f_ub, PK=PK, method="Brent",
          lower = 0, upper = 5)

#Version 3:
loss_v3 = function(bbeta, vv, f_lb, f_ub, PK) {
  (PK-preference_density(bbeta, vv, f_lb,f_ub))^2
}

ov3 = optim(0.2, loss_v3, vv=3, f_lb=f_lb, f_ub=f_ub, PK=PK, method="Brent",
          lower = 0, upper = 5)

print(c("beta_1:",ov1$par,
        "beta_2:",ov2$par,
        "beta_3:",ov3$par))


# Part D: Compare magnitude of elasticities:
#-------------------------------------------
# identification can't be achieved since changes in the distribution change 
# the elasticity \beta.


# Part E: Compute bounds on taxable income elasticity:
#-----------------------------------------------------
# Optimize lower bound

ssigma = 1

lower_bound_beta <- function(bbeta,sigma){
  D_left <- f_lb*(ub*(rho1/rho2)^bbeta - lb) 
  D_right <- f_ub*(ub - lb*(rho2/rho1)^bbeta) 
  Pk_hat <- sigma*max(D_left, D_right)
}

upper_bound_beta <- function(bbeta,sigma){
  D_left <- f_lb*(ub*(rho1/rho2)^bbeta - lb) 
  D_right <- f_ub*(ub - lb*(rho2/rho1)^bbeta) 
  Pk_hat <- sigma*min(D_left, D_right)
}

# Compute lower bound for beta
loss_lb <- function(bbeta, sigma, H_K){(H_K - lower_bound_beta(bbeta, sigma))^2}
beta_lb = optim(0.2, loss_lb, sigma=ssigma, H_K=PK, method="Brent",
            lower = 0, upper = 5)

# Compute upper bound for beta
loss_lb <- function(bbeta, sigma, H_K){(H_K - upper_bound_beta(bbeta, sigma))^2}
beta_ub = optim(0.2, loss_lb, sigma=ssigma, H_K=PK, method="Brent",
                lower = 0, upper = 5)

print(c(
"Beta Lower Bound",beta_lb[[1]],
"Beta Upper Bound",beta_ub[[1]]
))

# Part F: What condition on \beta would justify linearity assumption.
#--------------------------------------------------------------------

# Here we're just plugging our estimate of \beta_1 into the linear approximation
# of preferences distribution.

bbeta_1 <- ov1$par
P_K_vf<- 1/2*(f_lb*rho1^bbeta_1 
     + f_ub*rho2^bbeta_1)*(ub * rho2^(-bbeta_1) - lb * rho1^(-bbeta_1))



