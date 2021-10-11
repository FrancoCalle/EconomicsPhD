import numpy as np
import matplotlib.pyplot as plt

#Fix J_1

mu = 0.5
sigma = 1.2
I   = 10000
J1  = 4000
J0 = I-J1
np.random.seed(114234)
Y = np.random.normal(mu, sigma, I)

def randomAssignment(I, J1):
    arr = np.zeros(I)
    arr[:J1]  = 1
    np.random.shuffle(arr)
    return arr.astype(int)

# One shuffle:

def computeMeanDifferences(Y,D,I,J1):

    T = sum(np.multiply(D, Y))/J1 - sum(np.multiply(1-D, Y))/(I - J1)

    return T

def computeMeanDifferencesDistribution(Y,I,J1):

    D = randomAssignment(I, J1)
    T = computeMeanDifferences(Y,D,I,J1)

    return T

tDist = [computeMeanDifferencesDistribution(Y, I, J1) for i in range(1000)]
plt.hist(tDist, bins = 20)
plt.savefig('tStatisticDistribution4000.png')

tDist = [computeMeanDifferencesDistribution(Y, I, 1000) for i in range(1000)]
plt.hist(tDist, bins = 20)
plt.savefig('tStatisticDistribution1000.png')


tDist = [computeMeanDifferencesDistribution(Y, I, 200) for i in range(1000)]
plt.hist(tDist, bins = 20)
plt.savefig('tStatisticDistribution200.png')

tDist = [computeMeanDifferencesDistribution(Y, I, 100) for i in range(1000)]
plt.hist(tDist, bins = 20)
plt.savefig('tStatisticDistribution100.png')
