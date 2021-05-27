import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df_full = pd.read_stata('https://dl.dropboxusercontent.com/s/reqqx19nliefyw4/NLSY.dta')
df_extract = pd.read_stata('https://dl.dropboxusercontent.com/s/v7dt5q45eysau9y/nlsy79_extract.dta')
df_extract['male_raw'] = df_extract['male']
df_extract['male'] = (df_extract['male_raw'] == 'MALE')

df = df_full.copy()


## Data cleaning:
##---------------
df_educ = df_extract.groupby('id').apply(lambda df: df['yearseduc_age30'].mean())
df_educ.name = 'yearseduc_age30'
df = df.join(df_educ, on=['id'])

# Drop Missings:
df.dropna(
    subset=[
    'hrs_worked_pastyear',
    'yearly_inc',
    'nchld',
    'marstat',
    'yearseduc_age30',
    'family_net_wealth'
    ],
    inplace=True
)


# Define number of children, capped at 5
df['children_max5'] = df['nchld'].clip(upper=5)
# Define marital status indicators
df['marstat_2'] = (df['marstat'] == 2).astype(np.float64)
df['marstat_3'] = (df['marstat'] == 3).astype(np.float64)
# Define family net worth (in millions)
df['family_net_wealth_millions'] = df['family_net_wealth'] / 1e6

# Define people who worked as people who worked 100 or more hours (lower values
# are possibly misresponses or otherwise edge cases) AND who earned more than $500
# dollars (same rationale)

df['worked'] = (df['hrs_worked_pastyear'] > 100.) & (df['yearly_inc'] > 500.)
df['hourly_wage'] = (df['yearly_inc'] / (df['hrs_worked_pastyear'] + 1e-50)) *df['worked']
df['ln_hourly_wage'] = np.log(df['hourly_wage'])
df = df.loc[df['ln_hourly_wage']>-np.inf]

# People are identified by id:
df.sort_values('id').head()

# Check how many obs are by year:
df['year'].sort_values().value_counts()

# Check observations by id:
id_counts = (df['id'].value_counts()
            .reset_index()
            .rename(columns={'id': 'counts', 'index': 'id'})
            .query('counts == 16')
            .id
            )


nlsyBalancedDataFrame = df.loc[df['id'].isin(id_counts),:]

# Add year and individual dummies for fixed effects (Dropped first):
df_year_dummies = pd.get_dummies(
    nlsyBalancedDataFrame['year'].astype(int), drop_first=True, prefix='year', dtype=np.float64
)

df_individual_dummies = pd.get_dummies(
    nlsyBalancedDataFrame['id'].astype(int), drop_first=True, prefix='id', dtype=np.float64
)

year_dummy_columns = list(df_year_dummies.columns)[1:]
individual_dummy_columns = list(df_individual_dummies.columns)

# Sort data
nlsyBalancedDataFrame = pd.concat([nlsyBalancedDataFrame, df_year_dummies, df_individual_dummies], axis=1).sort_values(by = ['id', 'year'])

# Get lagged outcome and other covariates:
lagedDataFrame = nlsyBalancedDataFrame[['id',
                'year',
                'ln_hourly_wage',
                'age',
                'yearseduc_age30',
                'family_net_wealth_millions'
                ]].rename(columns={'ln_hourly_wage':'ln_hourly_wage_lag',
                'age':'age_lag',
                'yearseduc_age30':'yearseduc_age30_lag',
                'family_net_wealth_millions':'family_net_wealth_millions_lag'})

lagedDataFrame['year'] = (lagedDataFrame['year'] + 1)

#Obtain lagged variable:
nlsyBalancedDataFrame = nlsyBalancedDataFrame.merge(lagedDataFrame, how = 'left', on=['id','year'])

#Drop observations for first year
nlsyBalancedDataFrame = nlsyBalancedDataFrame.dropna(subset=['ln_hourly_wage_lag'])

nlsyBalancedDataFrame['constant'] = 1


## Construct matrices for estimation:
##----------------------------------
marketWageDeterminants = [
                        'constant', 'age', 'family_net_wealth_millions', #'yearseduc_age30',
                        ]

marketWageDeterminantsLag = [
                        'age_lag','family_net_wealth_millions_lag',
                        ]

laggedDependentVariable = ['ln_hourly_wage_lag']

dependentVariable = ['ln_hourly_wage']




# Get depvar into numpy for OLS:
# ------------------------------

yLog = nlsyBalancedDataFrame[dependentVariable].values
yLogLagged = nlsyBalancedDataFrame[laggedDependentVariable].values

# Define function to run ols regressions:
def olsRegression(XX, Y):

    bbeta = np.dot(np.linalg.inv(np.dot(XX.T,XX)), np.dot(XX.T,Y))

    return bbeta


# Part A: First Model Compute regression with fixed effects and lagged variable:
#-------------------------------------------------------------------------------

# Build XXs
XX = nlsyBalancedDataFrame[laggedDependentVariable + individual_dummy_columns  + marketWageDeterminants].values

β_a = olsRegression(XX, yLog)

predictionError1A = yLog - np.matmul(XX,β_a)

np.mean(predictionError1A**2)



# Part B: Second Model Compute regression with error autocorrelation (This looks Sexier!):
#-----------------------------------------------------------------------------------------

XX = nlsyBalancedDataFrame[laggedDependentVariable + individual_dummy_columns + marketWageDeterminants + marketWageDeterminantsLag].values

β_a = olsRegression(XX, yLog)

predictionError1B = yLog - np.matmul(XX,β_a)

np.mean(predictionError1B**2)



# Part C: Third Model Compute regression dependent variable first difference:
#----------------------------------------------------------------------------

XX = nlsyBalancedDataFrame[laggedDependentVariable + individual_dummy_columns + marketWageDeterminants + marketWageDeterminantsLag].values
YY = yLog - yLogLagged

β_a = olsRegression(XX, YY)

predictionError1C = yLog - np.matmul(XX,β_a)

np.mean(predictionError1C**2)


# Part D: Same as B but ρ = 1:
#----------------------------

XX = nlsyBalancedDataFrame[['age', 'family_net_wealth_millions']].values
XXLag = nlsyBalancedDataFrame[marketWageDeterminantsLag].values
DeltaXX = XX - XXLag

DeltaYY = yLog - yLogLagged

β_a = olsRegression(DeltaXX, DeltaYY)

predictionError1D = DeltaYY - np.matmul(DeltaXX,β_a)

np.mean(predictionError1D**2)

# ARMA section:
#--------------

from arma_mle import mlePanelARMA
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce

# Test Model E.A:
arma = mlePanelARMA()

# Arrange data for estimation:
#-----------------------------

y = (nlsyBalancedDataFrame
    .sort_values(by=['id','year'])
    .pivot_table(index = 'id',
                    columns='year',
                    values = dependentVariable)).values.T

x1 = (nlsyBalancedDataFrame
    .sort_values(by=['id','year'])
    .pivot_table(index = 'id',
                    columns='year',
                    values = marketWageDeterminants[1])).values.T

x2 = (nlsyBalancedDataFrame
    .sort_values(by=['id','year'])
    .pivot_table(index = 'id',
                    columns='year',
                    values = marketWageDeterminants[2])).values.T

X = [x1, x2]

# Part E: All of the above but with ϵ_it is MA(1):
#-------------------------------------------------
#-------------------------------------------------

# E.A: MLE including Fixed Effects
#---------------------------------

y_demean, xx_demean, YY_mean, XX_mean = arma.demeanVariables(y,X)

initial_parameters =  np.ones(5)*.5

fun = lambda parameters: arma.objectiveFunction(
y_demean,
xx_demean,
parameters,
arma.maximumLikelihood
)

# Optimize
result = arma.optimize(fun, initial_parameters)

#Retrieve Parameters
parameters_hat_A = result.x

#Predict model using estimates:
predictionErrorA = arma.predict_model(y_demean, xx_demean, parameters_hat_A)

y_hat = y_demean - predictionError + YY_mean
np.mean(predictionErrorA[1:,:].flatten()**2)




# E.B: Regression with error autocorrelation :
#---------------------------------------------

initial_parameters =  np.ones(8)*.5
initial_parameters[0] = parameters_hat_A[0]
initial_parameters[1] = parameters_hat_A[1]
initial_parameters[-1] = parameters_hat_A[-1]

fun = lambda parameters: arma.objectiveFunction(
y_demean,
xx_demean,
parameters,
arma.maximumLikelihoodModelB
)

result = arma.optimize(fun, initial_parameters)

parameters_hat_B = result.x

predictionErrorB = arma.predict_modelB(y_demean, xx_demean, parameters_hat_B)

np.mean(predictionErrorB[2:,:].flatten()**2)




# E.C: Regression with error autocorrelation :
#---------------------------------------------
initial_parameters = parameters_hat_B.copy()

fun = lambda parameters: arma.objectiveFunction(
y_demean,
xx_demean,
parameters,
arma.maximumLikelihoodModelC
)

result = arma.optimize(fun, initial_parameters)

parameters_hat_C = result.x

predictionErrorC = arma.predict_modelC(y_demean, xx_demean, parameters_hat_C)

np.mean(predictionErrorC[2:,:].flatten()**2)


# E.D: Regression with first differences for Y and X:
#---------------------------------------------------

initial_parameters = parameters_hat_C.copy()
initial_parameters = np.append(initial_parameters[:2],initial_parameters[5:])

fun = lambda parameters: arma.objectiveFunction(
    y,
    X,
    parameters,
    arma.maximumLikelihoodModelD
)

result = arma.optimize(fun, initial_parameters)

parameters_hat_D = result.x

predictionErrorD = arma.predict_modelD(y_demean, xx_demean, parameters_hat_D)

np.mean(predictionErrorD[2:,:].flatten()**2)


### Compute likelihood ratios (Confirm that we are using the same data for both models):

arma.computeLikelihoodRatio(predictionError1A, predictionError1B, 1, 1)




change = y[:,1].std()
rrho = 0.22
dynamic = []
parameter_update = rrho

for ii in range(10):
    parameter_update =  parameter_update * rrho
    change += parameter_update*change
    dynamic.append(change)



change = y[:,1].std()
rrho = 0.21
dynamic2 = []
parameter_update = rrho

for ii in range(10):
    parameter_update =  parameter_update * rrho
    change += parameter_update*change
    dynamic2.append(change)


change = y[:,1].std()
rrho = 1
dynamic3 = []
parameter_update = rrho

for ii in range(10):
    parameter_update =  parameter_update * rrho
    change += parameter_update*change
    dynamic3.append(change)

plt.plot(dynamic)
plt.plot(dynamic2)
plt.plot(dynamic3)
