from arma_mle import mlePanelARMA
import numpy as np
import matplotlib.pyplot as plt

arma = mlePanelARMA(N=200, T = 16, K = 3)

# Test Model E.A:
#----------------

y, X = arma.generateFakeData()

y_demean, xx_demean, YY_mean, XX_mean = arma.demeanVariables(y,X)

# Generate anonymous function:
fun = lambda parameters: arma.objectiveFunction(y_demean,
                                                xx_demean,
                                                parameters,
                                                arma.maximumLikelihood)

aalpha0 = [0, 0, 0]
bbeta0 = [0]
phi0 = [0]
initial_parameters = np.array(aalpha0+bbeta0+phi0+[1.3])


result = arma.optimize(fun, initial_parameters)

parameters_hat = result.x

f_i = arma.recoverFixedEffects(y, X, parameters_hat)

plt.scatter(parameters_hat, arma.true_parameters[:6], color = 'maroon')
plt.plot(np.linspace(0,2),np.linspace(0,2), color = 'gray', linewidth=2)
plt.show()
#Chef Kiss

plt.scatter(np.append(parameters_hat, f_i), arma.true_parameters,  facecolors='none', edgecolors='b', linewidth=0.3)
plt.scatter(parameters_hat, arma.true_parameters[:6], color = 'maroon', linewidth=2)
plt.plot(np.linspace(min(arma.true_parameters),max(arma.true_parameters)),
        np.linspace(min(arma.true_parameters),max(arma.true_parameters)),
        color = 'gray', linewidth=2)
plt.show()
#Chef Kiss

predictionError = arma.predict_model(y, X, parameters_hat)

y_hat = y - predictionError

#Trim first period and plot:
plt.hist(y_hat[1:,:].flatten(), bins=50, facecolor='none',edgecolor='maroon', label='Predicted')
plt.hist(y[1:,:].flatten(), bins=50, facecolor='none',edgecolor='gray', label='True')
plt.legend()


# Test Model E.B:
#----------------
y, X = arma.generateFakeDataModelB()
y_demean, xx_demean, YY_mean, XX_mean = arma.demeanVariables(y,X)

# Generate anonymous function:
fun = lambda parameters: arma.objectiveFunction(
                                y_demean,
                                xx_demean,
                                parameters,
                                arma.maximumLikelihoodModelB
                                )

aalpha0 = [0.3, 0.2, 1]
aalpha0Prime = [-0.1*0.3, -0.1*0.5, -0.1*2]
bbeta0 = [0.5]
phi0 = [-.3]
phi2 = [.3]

initial_parameters = np.array(
                        aalpha0+
                        aalpha0Prime+
                        bbeta0+
                        phi0+
                        phi2+
                        [1.3]
                        )

# Without any type of restriction:
result = arma.optimize(fun, initial_parameters)
parameters_hat = result.x

K = 3
aalpha = parameters_hat[0:K]
aalphaPrime = parameters_hat[K:2*K]
bbeta = parameters_hat[2*K]
pphi = parameters_hat[2*K+1]
pphi2 = parameters_hat[2*K+2]
ssigma = np.abs(parameters_hat[-1])

# Recover structurals (Works well, but I need
# initial values at least for alphas0)

rrho = bbeta
print(aalpha,'\n', aalphaPrime/rrho)
print([pphi+rrho],'\n',phi2/rrho)
#Chef Kiss


predictionError = arma.predict_modelB(y, X, parameters_hat)

y_hat = y - predictionError

#Trim first period and plot (Chef Kiss!):
plt.hist(y_hat[2:,:].flatten(), bins=50, facecolor='none',edgecolor='maroon', label='Predicted')
plt.hist(y[2:,:].flatten(), bins=50, facecolor='none',edgecolor='gray', label='True')
plt.legend()



# Test Model E.B (Restricted) Perhaps use moment conditions instead:
#-------------------------------------------------------------------
y, X = arma.generateFakeDataModelB()
y_demean, xx_demean, YY_mean, XX_mean = arma.demeanVariables(y,X)

# Generate anonymous function:
fun = lambda parameters: arma.objectiveFunction(
                                y_demean,
                                xx_demean,
                                parameters,
                                arma.maximumLikelihoodModelBRestricted
                                )

aalpha0 = [0.3, 0.5, 2]
rrho0 = [0.5]
pphi0 = [-.2]

initial_parameters = np.array(
                        aalpha0+
                        rrho0+
                        pphi0+
                        [1.3]
                        )

# Without any type of restriction:
result = arma.optimize(fun, initial_parameters)
parameters_hat = result.x

K=3
aalpha = initial_parameters[0:K]
bbeta = initial_parameters[K]
pphi = initial_parameters[K+1]
ssigma = np.abs(initial_parameters[-1])
