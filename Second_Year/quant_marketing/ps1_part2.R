library("bayesm")
library("Rcpp")
library('JWileymisc')
library('tidyverse')
library('dplyr')

repo_path = "C:/Users/franc/Dropbox/Franco Econ Phd/2 Second Year/Winter/Quantitative Marketing/ps1/"
data_path = "C:/Users/franc/Dropbox/Franco Econ Phd/2 Second Year/Winter/Quantitative Marketing/ps1/pbout_final.csv" 

# Pbout Data:
pbout_data <- read.table(data_path, sep=",", head=TRUE)
pbout_data$product <- rep(1:9, length(pbout_data$panelid)/9)



# Heterogeneity in preferences:
# ------------------------------

R = 5000
p = 9 
ncoef = 4
nlgt = length(unique(pbout_data$panelid)) 
list_hh <- unique(pbout_data$panelid)
Mcmc1 = list(R=R, keep=1)

fixed_effects = do.call(rbind, replicate(nrow(pbout_data)/p, diag(p), simplify=FALSE))

output <- c()

for (count in 1:5){
  
  # Demean Z which explains variation at individual level 
  nz = 1
  Z = matrix(rep(1, nlgt),ncol=nz) 
  Z = t(t(Z) - apply(Z,1,mean)) 
  
  ## simulate data
  lgtdata1 = NULL
  lgtdata2 = NULL
  lgtdata3 = NULL
  
  for (i in 1:nlgt) {
    
    pbout_ii <- pbout_data %>% filter(panelid == list_hh[i])
    
    X1 = cbind(pbout_ii$price, pbout_ii$feature, pbout_ii$display, pbout_ii$loyalty)
    
    X2 = cbind(pbout_ii$price, fixed_effects[pbout_data$panelid == list_hh[i],1:8])  
    
    X3 = cbind(pbout_ii$price, pbout_ii$loyalty-mean(pbout_ii$loyalty), fixed_effects[pbout_data$panelid == list_hh[i],1:8]) 
    
    y = matrix(pbout_ii$choice, ncol = p, byrow = TRUE)
    
    y = y[,1]
    
    lgtdata1[[i]] = list(y=y, X=X1)
    
    lgtdata2[[i]] = list(y=y, X=X2)
    
    lgtdata3[[i]] = list(y=y, X=X3)
    
  }
  
  # How to construct the priors?
  
  Prior1 = list(ncomp=count)  
  
  #Mcmc1 = list(R=R, keep=5)
  
  Data1 = list(p=p, lgtdata=lgtdata1, Z=Z) 
  Data2 = list(p=p, lgtdata=lgtdata2, Z=Z) 
  Data3 = list(p=p, lgtdata=lgtdata3, Z=Z) 
  
  ## fit model without sign constraints
  
  out1 = rhierMnlRwMixture(Data=Data1, Prior=Prior1, Mcmc=Mcmc1)
  out2 = rhierMnlRwMixture(Data=Data2, Prior=Prior1, Mcmc=Mcmc1)
  out3 = rhierMnlRwMixture(Data=Data3, Prior=Prior1, Mcmc=Mcmc1)
  
  ll_1 <- mean(out1$loglike)
  ll_2 <- mean(out2$loglike)
  ll_3 <- mean(out3$loglike)
  
  output <- c(output, ll_1, ll_2, ll_3)
  
}

rowMeans(out1$betadraw[1,,])  # get mean coefficients
mat= matrix(output, nrow= 5)
write.csv(as.data.frame(mat), file = file.path(repo_path,'loglik_fmm.csv'))

# Calculate the optimal pricing:
#-------------------------------

# Again obtain beta posterior distribution:
Data3 = list(p=p, lgtdata=lgtdata3, Z=Z)
out3 = rhierMnlRwMixture(Data=Data3, Prior=list(ncomp=5), Mcmc=Mcmc1)

# The best fitting model is the third model with five transition states
# Now optimize price vector that maximizes profits.

mc <- c()

for (i in 1:8){
  mc <- c(mc, 0.7 * mean(pbout_data$price[pbout_data$product == i]))
}

mc <- c(mc, 0)
product <- 1:9

mc_table <- cbind(product, mc)
pbout_data <- merge(pbout_data, mc_table,by="product")

# Obtain optimal pricing:
mean_price <- mc[1:8]/0.7 

coefs <- apply(out3$betadraw, c(1,2), mean)

dim(coefs)

hh_info <- cbind(list_hh, coefs)

colnames(hh_info) <- c("panelid", "beta1", "beta2", "beta3", 
                       "beta4", "beta5", "beta6", 
                       "beta7", "beta8", "beta9", 
                       "beta10")

pbout_hh <- merge(pbout_data, hh_info, by = "panelid") 



# Calculate the profits using the in-sample mean price:

hypo_price <- c(mean_price, 0)

product <- 1:9

temp_table <- cbind(product, hypo_price)

pbout_temp <- merge(pbout_hh, temp_table, by = "product") 

pbout_temp$hypo_u <-  (pbout_temp$beta1 * pbout_temp$hypo_price + 
                         (pbout_temp$loyalty-mean(pbout_temp$loyalty)) * pbout_temp$beta2 + 
                         fixed_effects[,1] * pbout_temp$beta3 + 
                         fixed_effects[,2] * pbout_temp$beta4 + 
                         fixed_effects[,3] * pbout_temp$beta5 + 
                         fixed_effects[,4] * pbout_temp$beta6 + 
                         fixed_effects[,5] * pbout_temp$beta7 + 
                         fixed_effects[,6] * pbout_temp$beta8 + 
                         fixed_effects[,7] * pbout_temp$beta9 + 
                         fixed_effects[,8] * pbout_temp$beta10)


pbout_temp <- pbout_temp %>% 
  group_by(panelid, date) %>%
  mutate(sum_u = sum(exp(hypo_u))) %>%
  mutate(hypo_share = exp(hypo_u) / (1 + sum_u))

pbout_temp$profits <- pbout_temp$hypo_share * (pbout_temp$hypo_price - pbout_temp$mc)

init_profits <- sum(pbout_temp$profits)




target_function <- function(parameters){
  
  hypo_price <- c(parameters, 0)
  product <- 1:9
  temp_table <- cbind(product, hypo_price)
  pbout_hhtemp <- merge(pbout_hh, temp_table, by = "product") 
  
  pbout_hhtemp$hypo_u <-  (pbout_hhtemp$beta1 * pbout_hhtemp$hypo_price + 
                             (pbout_hhtemp$loyalty-mean(pbout_hhtemp$loyalty)) * pbout_hhtemp$beta2 + 
                             fixed_effects[,1] * pbout_hhtemp$beta3 + 
                             fixed_effects[,2] * pbout_hhtemp$beta4 + 
                             fixed_effects[,3] * pbout_hhtemp$beta5 + 
                             fixed_effects[,4] * pbout_hhtemp$beta6 + 
                             fixed_effects[,5] * pbout_hhtemp$beta7 + 
                             fixed_effects[,6] * pbout_hhtemp$beta8 + 
                             fixed_effects[,7] * pbout_hhtemp$beta9 + 
                             fixed_effects[,8] * pbout_hhtemp$beta10)  
  
  pbout_hhtemp <- pbout_hhtemp %>% 
    group_by(panelid, date) %>%
    mutate(sum_u = sum(exp(hypo_u))) %>%
    mutate(hypo_share = exp(hypo_u) / (1 + sum_u))
  
  pbout_hhtemp$profits <- pbout_hhtemp$hypo_share * (pbout_hhtemp$hypo_price - pbout_hhtemp$mc)
  
  profits <- sum(pbout_hhtemp$profits)
  
  print(init_profits - profits)
  
  return(init_profits - profits)
  
}


result <- optim(par = mean_price, fn = target_function, method = "BFGS")

price_hat <- result$par

print("Optimized price vector:")

print(price_hat)


# Obtain posterior profits:

posterior_profits <- c()

for (i in 1:200){
  
  draw_i <- out3$betadraw[,,i]
  
  hh_info <- cbind(list_hh, draw_i)
  colnames(hh_info) <- c("panelid", "beta1", "beta2", "beta3", 
                         "beta4", "beta5", "beta6", 
                         "beta7", "beta8", "beta9", 
                         "beta10")
  pbout_hh <- merge(pbout_data, hh_info, by = "panelid") 
  
  optim_price <- c(price_hat, 0)
  product <- 1:9
  temp_table <- cbind(product, optim_price)
  pbout_temp <- merge(pbout_hh, temp_table, by = "product") 
  
  #pbout_temp$u <-  pbout_temp$beta1 * pbout_temp$optim_price + pbout_temp$feature * pbout_temp$beta2 + pbout_temp$display * pbout_temp$beta3 + pbout_temp$beta4*pbout_temp$loyalty
  pbout_temp$u <-  (pbout_temp$beta1 * pbout_temp$optim_price + 
                      (pbout_temp$loyalty-mean(pbout_temp$loyalty)) * pbout_temp$beta2 + 
                      fixed_effects[,1] * pbout_temp$beta3 + 
                      fixed_effects[,2] * pbout_temp$beta4 + 
                      fixed_effects[,3] * pbout_temp$beta5 + 
                      fixed_effects[,4] * pbout_temp$beta6 + 
                      fixed_effects[,5] * pbout_temp$beta7 + 
                      fixed_effects[,6] * pbout_temp$beta8 + 
                      fixed_effects[,7] * pbout_temp$beta9 + 
                      fixed_effects[,8] * pbout_temp$beta10)  
  
  pbout_temp <- pbout_temp %>% 
    group_by(panelid, date) %>%
    mutate(sum_u = sum(exp(u))) %>%
    mutate(hypo_share = exp(u) / (1 + sum_u))
  
  pbout_temp$profits <- pbout_temp$hypo_share * (pbout_temp$optim_price - pbout_temp$mc)
  
  diff_profits <- sum(pbout_temp$profits) - init_profits
  print(diff_profits)
  
  posterior_profits <- c(posterior_profits, diff_profits)
  
}


hist(posterior_profits, breaks=40)
png(file=file.path(repo_path,'profit_difference.png'),
    width=600, height=350)



