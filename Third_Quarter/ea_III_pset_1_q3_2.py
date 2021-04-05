import numpy as np
import matplotlib.pyplot as plt

N = 10000
S = 10000
sigma = 1

def generateData(sigma, N):

    cov = [[sigma, 0], [0, sigma]]
    U = np.random.multivariate_normal((0,0), cov, size=N, check_valid='warn', tol=1e-8)
    U0 = U[:,1]
    U1 = U[:,0]
    Y0 = 2 + U[:,1]
    Y1 = 5 + U[:,0]
    D  = np.random.binomial(1, .5, size= N)
    Y = Y0 + D * (Y1 - Y0)
    beta_true = np.mean(Y1 - Y0)

    return Y0, Y1, Y, D


def generateCovariates(N,D):

    D.shape[0]
    X = np.concatenate([np.ones(N).reshape(N,1), D.reshape(N,1)],1)

    return  X

def olsEstimation(Y,X):

    N = D.shape[0]
    beta_hat = np.dot(np.linalg.inv(np.dot(X.T, X)),np.dot(X.T,Y))
    sigma_hat = np.sqrt(np.mean((Y - np.dot(X,beta_hat.reshape(np.max(beta_hat.shape),1)))**2))
    se = np.sqrt(np.diag(sigma**2*np.linalg.inv(np.dot(X.T, X)))/N)

    return beta_hat, sigma_hat, se


# Part A:
Y0, Y1, Y, D = generateData(sigma, N)
X = generateCovariates(N,D)
beta_hat, sigma_hat, se = olsEstimation(Y.reshape(N,1),X)


# Part B: Bootstrap
X_sample_list = []
Y_sample_list = []
for s in range(S):
    index = np.random.choice(N, N) # This is automatically assuming uniform distribution 1/N
    Y_sample =  Y.reshape(N,1)[index,:]
    X_sample =  X[index,:]
    Y_sample_list.append(Y_sample)
    X_sample_list.append(X_sample)


beta_hat_list = []
for XX,YY in zip(X_sample_list,Y_sample_list):

    beta_hat, sigma_hat, se = olsEstimation(YY,XX)
    beta_hat_list.append(beta_hat.T)

beta_hat_list = np.concatenate(beta_hat_list,0)

plt.hist(beta_hat_list[:,0], bins = 60)
plt.show()
plt.savefig("q3_p2_Beta_1.png")


plt.hist(beta_hat_list[:,1], bins = 60)
plt.show()
plt.savefig("q3_p2_Beta_2.png")
