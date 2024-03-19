import numpy as np
import pandas as pd
import subprocess
import sys

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

install('pyblp')

import pyblp

pyblp.options.digits = 2
pyblp.options.verbose = False
pyblp.__version__

nevo_dataset = pd.read_csv(pyblp.data.NEVO_PRODUCTS_LOCATION)
nevo_dataset.columns




product_data = pd.read_csv('ps1_ex4.csv').rename(columns={'market': 'market_ids',
                                                          'p': 'prices',
                                                          'z1': 'demand_instruments0',
                                                          'z2': 'demand_instruments1',
                                                          'z3': 'demand_instruments2',
                                                          'z4': 'demand_instruments3',
                                                          'z5': 'demand_instruments4',
                                                          'z6': 'demand_instruments5',
                                                          'choice':'product_ids'})

X1_formulation = pyblp.Formulation('1 + prices + x', absorb='C(market_ids)')
X2_formulation = pyblp.Formulation('1 + prices + x')
product_formulations = (X1_formulation, X2_formulation)
product_formulations


mc_integration = pyblp.Integration('monte_carlo', size=50, specification_options={'seed': 0})
mc_integration


pr_integration = pyblp.Integration('product', size=5)
pr_integration


mc_problem = pyblp.Problem(product_formulations, product_data, integration=mc_integration)
mc_problem

pr_problem = pyblp.Problem(product_formulations, product_data, integration=pr_integration)
pr_problem

bfgs = pyblp.Optimization('bfgs', {'gtol': 1e-4})
bfgs

CovMatrix = np.ones((3, 3))
CovMatrix[-1,-1]=0

results1 = mc_problem.solve(sigma=CovMatrix, optimization=bfgs)
results1

results2 = pr_problem.solve(sigma=CovMatrix, optimization=bfgs)
results2

results3 = mc_problem.solve(sigma=CovMatrix, optimization=bfgs)
results3



# Demographics:

agent_data = pd.read_csv(pyblp.data.NEVO_AGENTS_LOCATION)
agent_data.head(40)

agent_formulation = pyblp.Formulation('0 + income + income_squared + age + child')
agent_formulation


# Problem:

nevo_problem = pyblp.Problem(
product_formulations,
product_data,
agent_formulation,
agent_data
)
nevo_problem


# Initial Conditions:

initial_sigma = np.diag([0.3302, 2.4526, 0.0163, 0.2441])
initial_pi = np.array([
                        [ 5.4819, 0, 0.2037, 0 ],
                        [15.8935, -1.2000, 0, 2.6342],
                        [-0.2506, 0, 0.0511, 0 ],
                        [ 1.2650, 0, -0.8091, 0 ]
                        ])
tighter_bfgs = pyblp.Optimization('bfgs', {'gtol': 1e-5})

nevo_results = nevo_problem.solve(
                                initial_sigma,
                                initial_pi,
                                optimization=tighter_bfgs,
                                method='1s' 
                )
nevo_results