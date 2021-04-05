import numpy as np
import matplotlib.pyplot as plt


def generateCovariates(N):

    X = np.concatenate([np.ones(N).reshape(N,1), np.random.normal(0, 1, N).reshape(N,1)],1)

    return  X

def generateData(sigma,beta_true,N, X):

    U = np.random.normal(0, sigma, N).reshape(N,1)
    Y = np.dot(X,beta_true.reshape(2,1)) + U

    return  U, Y

def olsEstimation(Y,X):

    beta_hat = np.dot(np.linalg.inv(np.dot(X.T, X)),np.dot(X.T,Y))
    sigma_hat = np.sqrt(np.mean((Y - np.dot(X,beta_hat.reshape(2,1)))**2))
    se = np.sqrt(np.diag(sigma**2*np.linalg.inv(np.dot(X.T, X)))/N)

    return beta_hat, sigma_hat, se


def parametric_bootstrap(X, sigma, beta_true, S, N):

    beta_hat_list = []
    for s in range(S):
        U, Y = generateData(sigma, beta_true, N, X)
        beta_hat, sigma_hat, se = olsEstimation(Y,X)
        beta_hat_list.append(beta_hat.T)

    beta_hat_list = np.concatenate(beta_hat_list,0)

    return beta_hat_list



# Question 3a:
beta_true = np.array([2,3])
mu = 0
sigma = 2
N = 10000
S = 10000

X = generateCovariates(N)
U, Y = generateData(sigma,beta_true,N, X)
beta_hat, sigma_hat, se = olsEstimation(Y,X)

# Question 3b:
beta_hat_list = parametric_bootstrap(X, sigma, beta_true, S, N)
plt.hist(beta_hat_list[:,0], bins = 60)
plt.show()
plt.savefig("q3_p1_Beta_2.png")
