# Quant Marketing Pset 1:

#install.packages("bayesm")

#install.packages("Rcpp")

library("bayesm")

library("Rcpp")

#git_path = "C:/Users/franc/OneDrive/Documents/GitHub/Chicagobooth/EconomicsPhD/Second_Year/quant_marketing"  

data_path = "C:/Users/franc/Dropbox/Franco Econ Phd/2 Second Year/Winter/Quantitative Marketing/ps1/pbout_final.csv" #df_name = "pbout_final.txt"


# Pbout Data:
pbout_data <- read.table(data_path, sep=",", head=TRUE)

p = 9

R <- 2000 # or if we want a short test, set R = 20

# Specification 1: Only Price:
# ----------------------------

y = matrix(pbout_data$choice, ncol = p, byrow = TRUE)[,1]

X = matrix(pbout_data$price, nrow = length(pbout_data$panelid), byrow = TRUE)  # here, k is 1, only price is the characteristics; p = 9 the number of choices

Data1 = list(y = y, X = X, p = p)

Mcmc1 = list(R=R, keep=1)

out = rmnlIndepMetrop(Data=Data1, Mcmc=Mcmc1) ####### This gives us the estimate for beta, the coefficient for x (price here) in the problem 

print(summary(out$betadraw), digits = 10)

posteriors <- out$betadraw

pos_ll1 <- logMargDenNR(out$loglike)


# Specification 2: Alternative specific dummy and price:
# ------------------------------------------------------


fixed_effects = do.call(rbind, replicate(nrow(pbout_data)/p, diag(p), simplify=FALSE))

X2= cbind(pbout_data$price, fixed_effects)  

Data2 = list(y = y, X = X2, p = p)

out2 = rmnlIndepMetrop(Data=Data2, Mcmc=Mcmc1)

print(summary(out2$betadraw), digits = 10)

pos_ll2 <- logMargDenNR(out2$loglike)


# Specification 3: Alternative specific dummy, price, promotions:
# ---------------------------------------------------------------

fixed_effects = do.call(rbind, replicate(nrow(pbout_data)/p, diag(p), simplify=FALSE))

X3= cbind(pbout_data$price, pbout_data$loyalty-mean(pbout_data$loyalty), fixed_effects[,1:8])  # , 

Data3 = list(y = y, X = X3, p = p)

out3 = rmnlIndepMetrop(Data=Data3, Mcmc=Mcmc1)

print(summary(out3$betadraw), digits = 10)

pos_ll3 <- logMargDenNR(out3$loglike)





