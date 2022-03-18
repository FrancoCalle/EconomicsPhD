library("bayesm")
library("Rcpp")
library('JWileymisc')
library('tidyverse')
library('dplyr')
library('fastDummies')

data_path = "C:/Users/franc/Dropbox/Franco Econ Phd/2 Second Year/Winter/Quantitative Marketing/ps1/pbout_final.csv" 

# Pbout Data:
pbout_data <- read.table(data_path, sep=",", head=TRUE)

pbout_data$product <- rep(1:9, length(pbout_data$panelid)/9)

R = 1000
p = 9 
nlgt = length(unique(pbout_data$panelid)) 
list_hh <- unique(pbout_data$panelid)
Mcmc1 = list(R=R, keep=1)


# Adding loyalty variable:
#-------------------------

pbout_brand <- dummy_cols(pbout_data, select_columns = "product")

intercept <- rep(1, nrow(pbout_data))

pbout_brand <- cbind(pbout_brand, intercept)


# 1. Only dummies:
#-----------------

y = matrix(pbout_brand$choice, ncol = p, byrow = TRUE)[,1]

X = cbind(pbout_brand$price, pbout_brand$loyalty)  

Data1 = list(y = y, X = X, p = p)

Mcmc1 = list(R=R, keep=1)

out0 = rmnlIndepMetrop(Data=Data1, Mcmc=Mcmc1) 

print(summary(out0$betadraw), digits = 10)

ll0 <- logMargDenNR(out0$loglike)

colMeans(out0$betadraw)




# 2. Brand intercept:
#--------------------

y = matrix(pbout_brand$choice, ncol = 9, byrow = TRUE)[,1]

X = cbind(pbout_brand$price, 
              pbout_brand$loyalty, 
              pbout_brand$product_1, 
              pbout_brand$product_2, 
              pbout_brand$product_3, 
              pbout_brand$product_4, 
              pbout_brand$product_5, 
              pbout_brand$product_6, 
              pbout_brand$product_7, 
              pbout_brand$product_8)  

Data1 = list(y = y, X = X, p = p)

Mcmc1 = list(R=R, keep=1)

out1 = rmnlIndepMetrop(Data=Data1, Mcmc=Mcmc1) 

print(summary(out1$betadraw), digits = 10)

ll1 <- logMargDenNR(out1$loglike)

colMeans(out1$betadraw)


# 3. Heterogeneity without branding:
#-----------------------------------

nlgt = length(unique(pbout_brand$panelid)) 

list_hh <- unique(pbout_brand$panelid)

nz = 1

Z = matrix(rep(1, nlgt),ncol=nz) 

Z = t(t(Z) - apply(Z,1,mean)) 

# Iterate accros costumers:

lgtdata = NULL

for (i in 1:nlgt) { 
  
  temp_data <- pbout_brand %>% filter(panelid == list_hh[i])
  
  X = cbind(temp_data$price, temp_data$loyalty)  
  
  y = matrix(temp_data$choice, ncol = p, byrow = TRUE)
  
  y = y[,1]
  
  lgtdata[[i]] = list(y=y, X=X)
}


Prior2 = list(ncomp=5) 

Mcmc2 = list(R=5000, keep=1)  

Data2 = list(p=p, lgtdata=lgtdata, Z=Z) 

out2 = rhierMnlRwMixture(Data=Data2, Prior=Prior2, Mcmc=Mcmc2)

ll2 <- mean(out2$loglike)

rowMeans(out2$betadraw[1,,])  



# 4. Heterogeneity and branding:
#-------------------------------

lgtdata2 = NULL

for (i in 1:nlgt) { 
  
  temp_data <- pbout_brand %>% filter(panelid == list_hh[i])
  
  X = cbind(temp_data$price, 
            temp_data$feature, 
            temp_data$display, 
            temp_data$loyalty, 
            temp_data$product_1, 
            temp_data$product_2, 
            temp_data$product_3, 
            temp_data$product_4, 
            temp_data$product_5, 
            temp_data$product_6, 
            temp_data$product_7, 
            temp_data$product_8)  
  
  y = matrix(temp_data$choice, ncol = p, byrow = TRUE)
  
  y = y[,1]
  
  lgtdata2[[i]] = list(y=y, X=X)
}


Prior3 = list(ncomp=5) 

Mcmc3 = list(R=R, keep=1)  

Data3 = list(p=p, lgtdata=lgtdata2, Z=Z) 

out3 = rhierMnlRwMixture(Data=Data3, Prior=Prior3, Mcmc=Mcmc3)

ll3 <- mean(out3$loglike)

rowMeans(out3$betadraw[1,,])


