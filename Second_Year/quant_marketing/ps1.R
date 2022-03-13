# Quant Marketing Pset 1:

install.packages("bayesm")
install.packages("Rcpp")

library("bayesm")
library("Rcpp")
library("")


git_path = "C:/Users/franc/OneDrive/Documents/GitHub/Chicagobooth/EconomicsPhD/Second_Year/quant_marketing"  

data_path = "C:/Users/franc/Dropbox/Franco Econ Phd/2 Second Year/Winter/Quantitative Marketing/ps1/pbout_final.csv"

df_name = "pbout_final.txt"


# Pbout Data:
pbout_data <- read.table(data_path, sep=",", head=TRUE)

#write.csv(pbout_data, 'C:/Users/franc/Dropbox/Franco Econ Phd/2 Second Year/Winter/Quantitative Marketing/ps1/pbout_2.csv')
#table(pbout_data["choice"]) # We have nine products ...



# llmnl(-999, data.matrix(pbout_data$choice), data.matrix(pbout_data$price))

y = pbout_data[,3]

X1 <- createX(p=9, na=1, Xa=pbout_data[,8:16], nd=NULL, Xd=NULL, base=1)


mydata <- list(y = y,
               X = X1,
               p = 9)

mymcmc <- list(R = 1000, nprint = 0)

out <- rmnlIndepMetrop(Data = mydata, Mcmc = mymcmc)

summary(out$betadraw)

plot(out$betadraw)


dim(out$betadraw)

median(out$betadraw[,1])
median(out$betadraw[,2])
median(out$betadraw[,3])
median(out$betadraw[,4])
median(out$betadraw[,5])
median(out$betadraw[,6])
median(out$betadraw[,7])
median(out$betadraw[,8])
median(out$betadraw[,9])



#------------------------------------------------------------------------------
# Testing model:
#------------------------------------------------------------------------------


data(margarine)
str(margarine)
marg <- merge(margarine$choicePrice, margarine$demos, by = "hhid")


y <- marg[,2]


X1 <- createX(p=10, na=1, Xa=marg[,3:12], nd=NULL, Xd=NULL, base=1)



out <- rmnlIndepMetrop(Data = list(y=y, X=X1, p=10), 
                       Mcmc = list(R=1e3, nprint=1e3))
















