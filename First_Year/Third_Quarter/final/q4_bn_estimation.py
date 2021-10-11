import numpy as np
import pandas as pd
import scipy
import statsmodels as sm
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

df_raw = pd.read_csv('https://dl.dropboxusercontent.com/s/9rp201wrk9m6qpb/data_final.csv',index_col=0)
columns = list(df_raw.columns)

df = df_raw[columns].copy()
df['D'] = df['D'].astype(np.float64)

k=11 # Treatment starts at 11
df_past_11 = df[df['t'] > 11].copy()
df_past_11.groupby('i')['D'].std().mean()

# Generate binary variable indicating if belongs to control or treatment group:
nObs = df[['D','i']].groupby('i').sum().shape[0]
ntreat = df[['D','i']].groupby('i').sum().query('D>0').shape[0]
shareTreat = ntreat/nObs
idTreatedList= df[['D','i']].groupby('i').sum().query('D>0').index.tolist() # Get list of treated people

# Add group tag:
df['Group: T'] = df['i'].isin(idTreatedList).astype(int)

# Define regression:
def olsRegression(XX, Y):

    bbeta = np.dot(np.linalg.inv(np.dot(XX.T,XX)), np.dot(XX.T,Y))

    return bbeta

# Get residual U:
df['U'] = df[['Y']].values - np.matmul(df[['X1','X2']].values,olsRegression(df[['X1','X2']].values,df[['Y']].values))

# Add Lagged data:
#-----------------

#First lag:
dfLag = df.copy().rename(columns={'Y': 'Y_Lag', 'U':'U_Lag', 'X1':'X1_Lag', 'X2':'X2_Lag', 'D':'D_Lag'})
dfLag['t'] = (dfLag['t'] + 1)

#Obtain lagged variable:
df = df.merge(dfLag[['Y_Lag','U_Lag','X1_Lag','X2_Lag','D_Lag','i','t']], how = 'left', on=['i','t'])

#Second lag:
dfLagLag = df.copy().rename(columns={'Y': 'Y_LagLag', 'U':'U_LagLag'})
dfLagLag['t'] = (dfLagLag['t'] + 2)

#Obtain lagged variable:
df = df.merge(dfLagLag[['Y_LagLag','U_LagLag','i','t']], how = 'left', on=['i','t'])

df_new = (df[['X1','X2','D','X1_Lag','X2_Lag','D_Lag','Y','Y_Lag','i','t','Group: T']]
            .copy()
            .dropna())

# Restrict data for observations at t ∈ [1, 30]:
df_new = (df[['X1','X2','D','X1_Lag','X2_Lag','D_Lag','Y','Y_Lag','i','t','Group: T']]
            .copy()
            .dropna()
            .query('`Group: T` == 1')
            .query('t >= 1')
            .query('t <= 30')
            )



zVariableList = ['X1','X2','D','X1_Lag','X2_Lag','D_Lag','Y_Lag']


def balestraNerloveEstimator(df, zVariableList):

    xi_mean_list=[]
    yi_mean_list=[]
    alpha_1_list=[]
    gamma_list = []
    id_list = df.i.unique()


    for ii in id_list:
        df_i = df.query('i =='+str(ii))
        xx_temp = df_i[zVariableList].values
        y_temp = df_i[['Y']].values
        gamma = olsRegression(xx_temp, y_temp)
        # Take Means:
        xx_temp_mean = xx_temp.mean(0).reshape(xx_temp.shape[1],1)
        y_temp_mean = y_temp.mean()
        xi_mean_list.append(xx_temp_mean)
        yi_mean_list.append(y_temp_mean)
        # Obtain alpha_1
        alpha_1 = y_temp_mean - np.matmul(gamma.T,xx_temp_mean)
        gamma = gamma.reshape(len(gamma),).tolist()
        alpha_1_list.append(alpha_1[0][0])
        gamma_list.append(gamma)

    gamma_hat = np.array(gamma_list).mean(0)
    alpha_1 = np.array(alpha_1_list).mean()

    #Compute sigma^2:
    Y_flat = df[['Y']].values
    X_flat = df[zVariableList].values
    uhat = (Y_flat - np.matmul(X_flat,
                                gamma_hat.reshape(len(gamma_hat),1)
                                ))
    sigma_residual = (1/df.shape[0]) * np.sum(uhat**2)

    #Compute sigma_u
    sigma_mu_residual = 0
    Nobs=len(xi_mean_list)
    for ii in range(len(xi_mean_list)):
        sigma_mu_residual += 1/Nobs * (yi_mean_list[ii] -
        np.matmul(gamma_hat,xi_mean_list[ii].reshape(xi_mean_list[ii].shape[0],1))
        - alpha_1)**2

    rho = sigma_mu_residual[0]/sigma_residual


    #varepsilon and eta:

    T = len(df.t.unique())
    varepsilon = sigma_residual*((1-rho)+T*rho)
    eta = sigma_residual*(1-rho)

    # Construct Omega matrix:
    C = np.eye(T)*eta
    C[0,0] = varepsilon
    I_n = np.eye(Nobs)
    OOmega = np.kron(I_n,C)

    # Transform x's and y's to do gls estimation:

    # Second Round Estimation:
    OOmega_inv = np.linalg.inv(OOmega)

    gammaGLS = np.dot(  np.linalg.inv(
                        np.dot(
                        np.dot(X_flat.T,OOmega_inv),
                        X_flat
                        )
                        ),
                        np.dot(
                        np.dot(X_flat.T,OOmega_inv),
                        Y_flat
                        )
                        )

    return gammaGLS

gammaGLS = balestraNerloveEstimator(df_new, zVariableList)
