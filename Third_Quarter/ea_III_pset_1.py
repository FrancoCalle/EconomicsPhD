import numpy as np
import matplotlib.pyplot as plt
from numpy.random import MT19937
from numpy.random import RandomState, SeedSequence
rs = RandomState(MT19937(SeedSequence(123456789)))

def generateData(sigma, rho, N):

    cov = [[sigma**2, rho * sigma], [rho * sigma, 1]]
    U = np.random.multivariate_normal((0,0), cov, size=N, check_valid='warn', tol=1e-8)
    U0 = U[:,1]
    U1 = U[:,0]
    D  = U1 > U0
    Y = D*U1 + (1-D)*U0

    return U0, U1, Y, D

def averageTreatmentEffect(sigma, rho, N, U1, U0, Y, D):

    ATE_estimated = np.mean(U1- U0)
    ATE_analytical = 0

    return ATE_estimated, ATE_analytical

def averageTreatmentOnTreated(sigma, rho, N, U1, U0, Y, D):

    ATT_estimated = -np.mean(U0[D == 1] - U1[D == 1])
    ATT_analytical = np.sqrt(2/np.pi)*np.sqrt(sigma**2 + 1 - 2*rho*sigma)

    return ATT_estimated, ATT_analytical


def averageTreatmentOnUntreated(sigma, rho, N, U1, U0, Y, D):

    ATU_estimated = np.mean(U1[D == 0] - U0[D == 0])
    ATU_analytical = -np.sqrt(2/np.pi)*np.sqrt(sigma**2 + 1 - 2*rho*sigma)

    return ATU_estimated, ATU_analytical

def beta_OLS(sigma, rho, N, U1, U0, Y, D):

    beta_estimated = np.mean(Y[D == 1]) - np.mean(Y[D == 0])
    beta_analytical = np.sqrt(2/np.pi)*(sigma**2 - 1)/np.sqrt(sigma**2 -2*rho*sigma +1)

    return beta_estimated, beta_analytical

# Model 1:

sigma = 2
rho = .5
N = 10000

U0, U1, Y, D = generateData(sigma, rho, N)
ATE_estimated1, ATE_analytical1 = averageTreatmentEffect(sigma, rho, N, U1, U0, Y, D)
ATT_estimated1, ATT_analytical1 = averageTreatmentOnTreated(sigma, rho, N, U1, U0, Y, D)
ATU_estimated1, ATU_analytical1 = averageTreatmentOnUntreated(sigma, rho, N, U1, U0, Y, D)
beta_estimated1, beta_analytical1 = beta_OLS(sigma, rho, N, U1, U0, Y, D)

print("\nModel 1: Sigma 2 and rho 0.5")
print("ATE Estimated:",ATE_estimated1,"\nATE Analytical:", ATE_analytical1)
print("ATT Estimated:",ATT_estimated1,"\nATT Analytical:", ATT_analytical1)
print("ATU Estimated:",ATU_estimated1,"\nATU Analytical:", ATU_analytical1)
print("Beta Estimated:",beta_estimated1,"\nBeta Analytical:", beta_analytical1)

# Model 2:

sigma = 2
rho = 0

U0, U1, Y, D = generateData(sigma, rho, N)
ATE_estimated2, ATE_analytical2 = averageTreatmentEffect(sigma, rho, N, U1, U0, Y, D)
ATT_estimated2, ATT_analytical2 = averageTreatmentOnTreated(sigma, rho, N, U1, U0, Y, D)
ATU_estimated2, ATU_analytical2 = averageTreatmentOnUntreated(sigma, rho, N, U1, U0, Y, D)
beta_estimated2, beta_analytical2 = beta_OLS(sigma, rho, N, U1, U0, Y, D)

print("\nModel 2: Sigma 2 and rho 0")
print("ATE Estimated:",ATE_estimated2,"\nATE Analytical:", ATE_analytical2)
print("ATT Estimated:",ATT_estimated2,"\nATT Analytical:", ATT_analytical2)
print("ATU Estimated:",ATU_estimated2,"\nATU Analytical:", ATU_analytical2)
print("Beta Estimated:",beta_estimated2,"\nBeta Analytical:", beta_analytical2)

# Model 3:

sigma = 2
rho = -.5

U0, U1, Y, D = generateData(sigma, rho, N)
ATE_estimated3, ATE_analytical3 = averageTreatmentEffect(sigma, rho, N, U1, U0, Y, D)
ATT_estimated3, ATT_analytical3 = averageTreatmentOnTreated(sigma, rho, N, U1, U0, Y, D)
ATU_estimated3, ATU_analytical3 = averageTreatmentOnUntreated(sigma, rho, N, U1, U0, Y, D)
beta_estimated3, beta_analytical3 = beta_OLS(sigma, rho, N, U1, U0, Y, D)

print("\nModel 3: Sigma 2 and rho -0.5")
print("ATE Estimated:",ATE_estimated3,"\nATE Analytical:", ATE_analytical3)
print("ATT Estimated:",ATT_estimated3,"\nATT Analytical:", ATT_analytical3)
print("ATU Estimated:",ATU_estimated3,"\nATU Analytical:", ATU_analytical3)
print("Beta Estimated:",beta_estimated3,"\nBeta Analytical:", beta_analytical3)


# Variation in Sigma:

rho = .5
sigma_list = np.linspace(1,3)
ATE_estimated_list = []
ATT_estimated_list = []
ATU_estimated_list = []
beta_estimated_list = []

for ssigma in sigma_list:

    U0, U1, Y, D = generateData(ssigma, rho, N)
    ATE_estimated, _ = averageTreatmentEffect(sigma, rho, N, U1, U0, Y, D)
    ATE_estimated_list.append(ATE_estimated)

    ATT_estimated, _ = averageTreatmentOnTreated(sigma, rho, N, U1, U0, Y, D)
    ATT_estimated_list.append(ATT_estimated)

    ATU_estimated, _ = averageTreatmentOnUntreated(sigma, rho, N, U1, U0, Y, D)
    ATU_estimated_list.append(ATU_estimated)

    beta_estimated, _ = beta_OLS(sigma, rho, N, U1, U0, Y, D)
    beta_estimated_list.append(beta_estimated)



plt.plot(sigma_list,ATT_estimated_list)
plt.savefig("ATT_sigma_test.png")

plt.plot(sigma_list,ATU_estimated_list)
plt.savefig("ATU_sigma_test.png")
