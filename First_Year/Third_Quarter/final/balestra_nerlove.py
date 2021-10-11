import numpy as np
import pandas as pd
import scipy
import statsmodels as sm
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from ar_mle import mlePanelAR

# Test Model E.A:
ar = mlePanelAR()

y, X = ar.generateFakeDataModelB()

def olsRegression(XX, Y):

    bbeta = np.dot(np.linalg.inv(np.dot(XX.T,XX)), np.dot(XX.T,Y))

    return bbeta


y_lag = y[1:,:]
XX_lag = [x[1:,:] for x in X]

y = np.delete(y,(-1), axis=0)
XX = [np.delete(x,(-1), axis=0) for x in X]

# Estimation procedure:
Nobs = y.shape[1]
T = y.shape[0]
gamma_list = []
alpha_1_list = []

yi_mean_list = []
xi_mean_list = []

for ii in range(y.shape[1]):
    xx_temp = np.array([
                        XX[0][:,ii],
                        XX[1][:,ii],
                        XX[2][:,ii],
                        XX_lag[0][:,ii],
                        XX_lag[1][:,ii],
                        XX_lag[2][:,ii],
                        y_lag[:,ii]]
                        )
    y_temp = y[:,ii].reshape(y.shape[0],1)
    gamma = olsRegression(xx_temp.T, y_temp)
    # Take Means:
    xx_temp_mean = xx_temp.mean(1).reshape(xx_temp.shape[0],1)
    y_temp_mean = y[:,ii].mean()
    xi_mean_list.append(xx_temp_mean)
    yi_mean_list.append(y_temp_mean)
    # Obtain alpha_1
    alpha_1 = y[:,ii].mean() - np.matmul(gamma.T,xx_temp_mean)
    gamma = gamma.reshape(len(gamma),).tolist()
    alpha_1_list.append(alpha_1[0][0])
    gamma_list.append(gamma)

gamma_hat = np.array(gamma_list).mean(0)
alpha_1 = np.array(alpha_1_list).mean()

#compute sigma_u and sigma
tt = 1
sigma_mu_residual = 0
sigma_residual = 0
for ii in range(y.shape[1]):
    sigma_mu_residual += 1/Nobs * (yi_mean_list[ii] -
    np.matmul(gamma_hat,xi_mean_list[ii].reshape(xi_mean_list[ii].shape[0],1))
    - alpha_1)**2
    for tt in range(y.shape[0]):
        sigma_residual += (
        1/(Nobs*T) * (y[tt,ii] -
        gamma_hat[0]*XX[0][tt,ii] -
        gamma_hat[1]*XX[1][tt,ii] -
        gamma_hat[2]*XX[2][tt,ii] -
        gamma_hat[3]*XX_lag[0][tt,ii] -
        gamma_hat[4]*XX_lag[1][tt,ii] -
        gamma_hat[5]*XX_lag[2][tt,ii] -
        gamma_hat[6]*y_lag[tt,ii] -
        alpha_1)**2
        )

rho = sigma_mu_residual[0]/sigma_residual

#varepsilon and eta:

varepsilon = sigma_residual*((1-rho)+T*rho)
eta = sigma_residual*(1-rho)

# Construct Omega matrix:
C = np.eye(T)*eta
C[0,0] = varepsilon
I_n = np.eye(Nobs)

OOmega = np.kron(I_n,C)

# Transform x's and y's to do gls estimation:
y_flat = y.flatten(order='F')
z_flat = np.array([x.flatten(order='F') for x in XX] +
                    [x.flatten(order='F') for x in XX_lag] +
                    [y_lag.flatten(order='F')]).T

# Second Round Estimation:
OOmega_inv = np.linalg.inv(OOmega)

gammaGLS = np.dot(  np.linalg.inv(
                    np.dot(
                    np.dot(z_flat.T,OOmega_inv),
                    z_flat
                    )
                    ),
                    np.dot(
                    np.dot(z_flat.T,OOmega_inv),
                    y_flat
                    )
                    )

gammaOLS = np.dot(  np.linalg.inv(
                    np.dot(z_flat.T,
                    z_flat
                    )
                    ),
                    np.dot(z_flat.T,
                    y_flat
                    )
                    )


print('GLS:', gammaGLS)
print('OLS:', gammaOLS)
print('HAT:', gamma_hat)


gammaGLS[-1]

gammaOLS[-1]
