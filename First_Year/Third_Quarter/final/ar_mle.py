import numpy as np
from scipy.optimize import minimize
from functools import reduce

class mlePanelAR:


    def __init__(self, N=20, T = 200, K = 3):

        # Global string variables for paths and relative paths:
        self.N = N
        self.T = T
        self.K = K
        self.seed = 1034248

        #Numerical optimization parameters:
        self.maxiter = 10000
        self.tolerance = 1e-15
        self.display = True
        self.method =  'nelder-mead'

        #Set True parameters for fake data:
        self.aalpha = [0.3, 0.5, .2]
        self.bbeta = 0.8
        self.rrho = 0.8
        self.sigma = 1.5
        self.mu = 0

        #Parameters for fixed effects:
        self.mu_fe = 0
        self.sigma_fe =2

        np.random.seed(self.seed)
        self.fixedEffects = np.random.normal(self.mu_fe,
                        self.sigma_fe,
                        self.N
                        )

        self.true_parameters = np.append(
        np.array(
        self.aalpha +
        [self.bbeta] +
        [self.sigma]
        ),
        self.fixedEffects
        )

        #Initial parameters:
        self.initial_parameters = np.ones(6)*.5


    def demeanVariables(self, y, X):

        YY = y - y.mean(0)
        XX = [X[kk] - X[kk].mean(0) for kk in range(len(X))]

        YY_mean = y.mean(0)
        XX_mean = [X[kk].mean(0) for kk in range(len(X))]

        return YY, XX, YY_mean, XX_mean


    def generateFakeData(self):

        #Generate XXs covariates and iid errors for each it:
        np.random.seed(self.seed)
        X = [np.random.normal(0,self.sigma,(self.T,self.N)) for kk in range(self.K)]
        epsilon_it = np.random.normal(self.mu,self.sigma,(self.T,self.N))
        epsilon_it[0,:] = 0 # First error term equals zero (otherwise identification unfeasible)

        y = np.zeros((self.T,self.N))

        # Compute yy for each time
        for ii in range(self.N):
            f_i = self.fixedEffects[ii]
            for tt in range(self.T-1):
                y[tt+1, ii] = (
                    f_i*0 +
                    self.bbeta * y[tt,ii] +
                    self.aalpha[0] * X[0][tt+1,ii] +
                    self.aalpha[1] * X[1][tt+1,ii] +
                    self.aalpha[2] * X[2][tt+1,ii] +
                    epsilon_it[tt+1, ii]
                    )

        y = y[1:,:]
        X[0] = X[0][1:,:]
        X[1] = X[1][1:,:]
        X[2] = X[2][1:,:]

        return y, X


    def generateFakeDataModelB(self):

        #Generate XXs covariates and iid errors for each it:
        np.random.seed(self.seed)
        X = [np.random.normal(0,self.sigma,(self.T,self.N)) for kk in range(self.K)]
        epsilon_it = np.random.normal(self.mu,self.sigma,(self.T,self.N))
        epsilon_it[0:2,:] = 0 # First error term equals zero (otherwise identification unfeasible)

        y = np.zeros((self.T,self.N))

        # Compute yy for each time
        for ii in range(self.N):
            f_i = self.fixedEffects[ii]
            for tt in range(self.T-2):
                y[tt+2, ii] = (
                    f_i*0 +
                    self.rrho * y[tt+1,ii]
                    + self.aalpha[0] * X[0][tt+2,ii]
                    + self.aalpha[1] * X[1][tt+2,ii]
                    + self.aalpha[2] * X[2][tt+2,ii]
                    - self.rrho * self.aalpha[0] * X[0][tt+1,ii]
                    - self.rrho * self.aalpha[1] * X[1][tt+1,ii]
                    - self.rrho * self.aalpha[2] * X[2][tt+1,ii]
                    + epsilon_it[tt+2, ii]
                    )

        y = y[2:,:]
        X[0] = X[0][2:,:]
        X[1] = X[1][2:,:]
        X[2] = X[2][2:,:]

        return y, X


    def normal_density(self, epsilon, sigma):

        psi = 1/(sigma * np.sqrt(2*np.pi))*np.exp(-1/2 * (epsilon/sigma)**2)

        return psi


    def maximumLikelihood(self, YY, XX, parameters):

        T = YY.shape[0]
        N = YY.shape[1]
        K = len(XX)
        aalpha = parameters[0:K]
        bbeta = parameters[K]
        ssigma = np.abs(parameters[-1])

        logL  = 0
        for ii in range(N):
            for tt in range(T-1):
                Z = bbeta * YY[tt,ii]  # First AR part
                for k in range(K): Z += aalpha[k] * XX[k][tt+1,ii] # Add covariates
                epsilon = YY[tt+1,ii] - Z

                psi = self.normal_density(epsilon, ssigma)
                logL += np.log(psi)

        return logL


    def maximumLikelihoodModelB(self, YY, XX, parameters):

        T = YY.shape[0]
        N = YY.shape[1]
        K = len(XX)
        aalpha = parameters[0:K]
        aalphaPrime = parameters[K:2*K]
        bbeta = parameters[2*K]
        ssigma = np.abs(parameters[-1])

        logL  = 0
        for ii in range(N):
            for tt in range(T-1):
                Z = bbeta * YY[tt,ii]  # First AR part
                for k in range(K): Z += aalpha[k] * XX[k][tt+1,ii] # Add time t covariates
                for k in range(K): Z -= aalphaPrime[k] * XX[k][tt,ii] # Add time t-1 covariates
                epsilon = YY[tt+1,ii] - Z
                psi = self.normal_density(epsilon, ssigma)
                logL += np.log(psi)

        return logL


    def maximumLikelihoodModelC(self, YY, XX, parameters):

        T = YY.shape[0]
        N = YY.shape[1]
        K = len(XX)
        aalpha = parameters[0:K]
        aalphaPrime = parameters[K:2*K]
        bbeta = parameters[2*K]
        ssigma = np.abs(parameters[-1])

        logL  = 0
        for ii in range(N):
            for tt in range(T-1):
                Z = bbeta * YY[tt,ii]  # First AR part
                for k in range(K): Z += aalpha[k] * XX[k][tt+1,ii] # Add time t covariates
                for k in range(K): Z -= aalphaPrime[k] * XX[k][tt,ii] # Add time t-1 covariates
                epsilon = (YY[tt+1,ii] - YY[tt,ii]) - Z
                psi = self.normal_density(epsilon, ssigma)
                logL += np.log(psi)

        return logL

    def maximumLikelihoodModelD(self, YY, XX, parameters):

        T = YY.shape[0]
        N = YY.shape[1]
        K = len(XX)
        aalpha = parameters[0:K]
        ssigma = np.abs(parameters[-1])

        logL  = 0
        for ii in range(N):
            for tt in range(T-1):
                Z = 0 * YY[tt,ii]  # First AR part
                for k in range(K): Z += aalpha[k] * (XX[k][tt+1,ii] - XX[k][tt,ii]) # Add time t covariates
                epsilon = (YY[tt+1,ii] - YY[tt,ii]) - Z
                psi = self.normal_density(epsilon, ssigma)
                logL += np.log(psi)

        return logL


    def objectiveFunction(self, YY, XX, parameters, maximumLikelihood):

        logL = maximumLikelihood(YY,XX, parameters)

        return -logL

    def optimize(self, anonymousFunction, initial_parameters, optimizationMethod=None):

        if optimizationMethod is None: optimizationMethod = self.method

        res = minimize(
                    anonymousFunction,
                    initial_parameters,
                    method=optimizationMethod,
                    options={'xatol': self.tolerance, 'disp': True, 'maxiter': 10000}
                    )

        return res

    def recoverFixedEffects(self, YY, XX, parameters):

        T = YY.shape[0]
        N = YY.shape[1]
        K = len(XX)

        bbeta = parameters[K]

        YY_mean = YY.mean(0)
        XX_mean = [XX[kk].mean(0) for kk in range(len(XX))]

        f_i = (1 - bbeta) * YY_mean


        return f_i


    def predict_model(self, YY, XX, parameters):

        T = YY.shape[0]
        N = YY.shape[1]
        K = len(XX)
        aalpha = parameters[0:K]
        bbeta = parameters[K]
        ssigma = np.abs(parameters[-1])

        predictionError = np.zeros((T, N))
        # First recover e_it
        for ii in range(N):
            for tt in range(T-1):
                Z = bbeta * YY[tt,ii]  # First AR part
                for k in range(K): Z += aalpha[k] * XX[k][tt+1,ii] # Add covariates
                epsilon = YY[tt+1,ii] - Z
                predictionError[tt+1,ii] = epsilon

        predictionError = predictionError[1:,:]

        return predictionError


    def predict_modelB(self, YY, XX, parameters):

        T = YY.shape[0]
        N = YY.shape[1]
        K = len(XX)
        aalpha = parameters[0:K]
        aalphaPrime = parameters[K:2*K]
        bbeta = parameters[2*K]
        ssigma = np.abs(parameters[-1])
        predictionError = np.zeros((T, N))

        for ii in range(N):
            for tt in range(T-1):
                Z = bbeta * YY[tt,ii]  # First AR part
                for k in range(K): Z += aalpha[k] * XX[k][tt+1,ii] # Add time t covariates
                for k in range(K): Z -= aalphaPrime[k] * XX[k][tt,ii] # Add time t-1 covariates
                epsilon = YY[tt+1,ii] - Z
                predictionError[tt+1,ii] = epsilon

        predictionError = predictionError[1:,:]

        return predictionError


    def predict_modelC(self, YY, XX, parameters):

        T = YY.shape[0]
        N = YY.shape[1]
        K = len(XX)
        aalpha = parameters[0:K]
        aalphaPrime = parameters[K:2*K]
        bbeta = parameters[2*K]
        ssigma = np.abs(parameters[-1])
        predictionError = np.zeros((T, N))

        for ii in range(N):
            for tt in range(T-1):
                Z = bbeta * YY[tt,ii]  # First AR part
                for k in range(K): Z += aalpha[k] * XX[k][tt+1,ii] # Add time t covariates
                for k in range(K): Z -= aalphaPrime[k] * XX[k][tt,ii] # Add time t-1 covariates
                epsilon = (YY[tt+1,ii] - YY[tt,ii]) - Z
                predictionError[tt+1,ii] = epsilon

        predictionError = predictionError[1:,:]

        return predictionError


    def predict_modelD(self, YY, XX, parameters):

        T = YY.shape[0]
        N = YY.shape[1]
        K = len(XX)
        aalpha = parameters[0:K]
        ssigma = np.abs(parameters[-1])
        predictionError = np.zeros((T, N))

        for ii in range(N):
            for tt in range(T-1):
                Z = 0 * YY[tt,ii]  # First AR part
                for k in range(K): Z += aalpha[k] * (XX[k][tt+1,ii]-XX[k][tt,ii]) # Add time t covariates
                epsilon = (YY[tt+1,ii] - YY[tt,ii]) - Z
                predictionError[tt+1,ii] = epsilon

        predictionError = predictionError[1:,:]

        return predictionError


    def computeLikelihoodRatio(self, error1, error2, sigma1, sigma2):

        pr1 = list(map(lambda x: np.log(self.normal_density(x,sigma1)), error1))
        pr2 = list(map(lambda x: np.log(self.normal_density(x,sigma2)), error2))

        L1 = reduce((lambda x, y: x + y), pr1)
        L2 = reduce((lambda x, y: x + y), pr2)

        LR = -2*(L2 - L1)

        return LR
