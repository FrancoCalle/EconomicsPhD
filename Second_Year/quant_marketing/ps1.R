# Quant Marketing Pset 1:

#install.packages("bayesm")
#install.packages("Rcpp")
#install.packages('JWileymisc')

library("bayesm")
library("Rcpp")
library('JWileymisc')

data_path = "C:/Users/franc/Dropbox/Franco Econ Phd/2 Second Year/Winter/Quantitative Marketing/ps1/pbout_final.csv" #df_name = "pbout_final.txt"


# Pbout Data:
pbout_data <- read.table(data_path, sep=",", head=TRUE)

p = 9

R <- 3000 # or if we want a short test, set R = 20

Mcmc1 = list(R=R, keep=1)



# Summary Statistics:
#-------------------

# First we create a summary statistics table
price <- matrix(pbout_data$price, nrow = length(pbout_data$panelid)/9, byrow = TRUE)
display <- matrix(pbout_data$display, nrow = length(pbout_data$panelid)/9, byrow = TRUE)
loyalty <- matrix(pbout_data$loyalty, nrow = length(pbout_data$panelid)/9, byrow = TRUE)

mean_price <- colMeans(price)[1:8]
mean_display <- colMeans(display)[1:8]
mean_loyalty <- colMeans(loyalty)[1:8]
choice_prob <- matrix(pbout_data$choice, nrow = length(pbout_data$panelid)/9, byrow = TRUE)[,1]

mean_choice_prob <- c(sum(choice_prob==1)/length(choice_prob), 
                      sum(choice_prob==2)/length(choice_prob),
                      sum(choice_prob==3)/length(choice_prob),
                      sum(choice_prob==4)/length(choice_prob),
                      sum(choice_prob==5)/length(choice_prob),
                      sum(choice_prob==6)/length(choice_prob),
                      sum(choice_prob==7)/length(choice_prob),
                      sum(choice_prob==8)/length(choice_prob))

mat = matrix(, nrow = 4, ncol = 8)
mat[1,] <- mean_price
mat[2,] <- mean_display*100
mat[3,] <- mean_loyalty*100
mat[4,] <- mean_choice_prob*100

print(round(mat,2))





# Specification 1: Only Price:
# ----------------------------

y = matrix(pbout_data$choice, ncol = p, byrow = TRUE)[,1]

X = matrix(pbout_data$price, nrow = length(pbout_data$panelid), byrow = TRUE)  # here, k is 1, only price is the characteristics; p = 9 the number of choices

Data1 = list(y = y, X = X, p = p)

out = rmnlIndepMetrop(Data=Data1, Mcmc=Mcmc1) ####### This gives us the estimate for beta, the coefficient for x (price here) in the problem 

print(summary(out$betadraw), digits = 10)

posteriors <- out$betadraw

quantile(out$betadraw, .05)
quantile(out$betadraw, .5)
quantile(out$betadraw, .95)

# Specification 2: Alternative specific dummy and price:
# ------------------------------------------------------


fixed_effects = do.call(rbind, replicate(nrow(pbout_data)/p, diag(p), simplify=FALSE))

X2= cbind(pbout_data$price, fixed_effects[,1:8])  

Data2 = list(y = y, X = X2, p = p)

out2 = rmnlIndepMetrop(Data=Data2, Mcmc=Mcmc1)

print(summary(out2$betadraw), digits = 10)

quantile(out2$betadraw[,1], .05)
quantile(out2$betadraw[,1], .5)
quantile(out2$betadraw[,1], .95)


# Specification 3: Alternative specific dummy, price, promotions:
# ---------------------------------------------------------------

fixed_effects = do.call(rbind, replicate(nrow(pbout_data)/p, diag(p), simplify=FALSE))

X3= cbind(pbout_data$price, pbout_data$loyalty-mean(pbout_data$loyalty), fixed_effects[,1:8])  # , 

Data3 = list(y = y, X = X3, p = p)

out3 = rmnlIndepMetrop(Data=Data3, Mcmc=Mcmc1)

print(summary(out3$betadraw), digits = 10)

quantile(out3$betadraw[,1], .05)
quantile(out3$betadraw[,1], .5)
quantile(out3$betadraw[,1], .95)


# Compute credibility interval:

u_q <- quantile(out$betadraw, 0.95)
l_q <- quantile(out$betadraw, 0.05)

u_q <- quantile(out2$betadraw[,1], 0.95)
l_q <- quantile(out2$betadraw[,1], 0.05)

u_q <- quantile(out3$betadraw[,1], 0.95)
l_q <- quantile(out3$betadraw[,1], 0.05)

# Compute posterior log likelihood:
pos_ll1 <- logMargDenNR(out$loglike[out$loglike < quantile(out$loglike, 0.95) & out$loglike > quantile(out$loglike, 0.05)])
pos_ll2 <- logMargDenNR(out2$loglike[out2$loglike < quantile(out2$loglike, 0.95) & out2$loglike > quantile(out2$loglike, 0.05)])
pos_ll3 <- logMargDenNR(out3$loglike[out3$loglike < quantile(out3$loglike, 0.95) & out3$loglike > quantile(out3$loglike, 0.05)])







