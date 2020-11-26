import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

class dynamicProgramingModule:


    def __init__(self, bbeta, z):

        # Global string variables for paths and relative paths:
        self.z = z
        self.bbeta =  bbeta
        self.ddelta = 0.1
        self.aalpha = 0.3
        self.gridSize = 500
        self.tol = 1e-5
        self.error = 500

    def BellmanEquation(self, k, y, V0):

        V = np.log(self.z*k**self.aalpha + (1-self.ddelta)*k - y) + self.bbeta*V0

        return V

    def smartInitialGuessV0(self, y):

        V_old = np.log(self.z*y**self.aalpha - self.ddelta*y)/(1-self.bbeta)

        return V_old

    def steadyStateVariables(self):

        #Actual Steady State of the model:
        steadyStateCapital = ((self.bbeta * self.z * self.aalpha)/(1- self.bbeta + self.ddelta*self.bbeta))**(1/(1-self.aalpha))
        steadyStateConsumption = self.z * steadyStateCapital**self.aalpha  - self.ddelta*steadyStateCapital

        return steadyStateCapital, steadyStateConsumption


    def find_nearest(self, a, a0):
        "Element in nd array `a` closest to the scalar value `a0`"
        idx = np.abs(a - a0).argmin()
        return a.flat[idx]


    def valueFunctionIteration(self, V_old, y):

        fig, ax = plt.subplots()
        ax.plot(y, V_old, color=plt.cm.jet(0), lw=1, alpha=0.6)

        error = self.error
        count = 0
        policyFunction = np.zeros(V_old.shape)

        while error>self.tol:

            V_new = np.ones(y.shape[0])

            for i in range(y.shape[0]):
                k = y[i]
                upperBoundY = (self.z*k**self.aalpha + (1-self.ddelta)*k)
                objective = [self.BellmanEquation(k,y[j],V_old[j]) if y[j] <= upperBoundY else -99999 for j in range(y.shape[0])]
                V_new[i] = np.max(objective)
                policyFunction[i] = y[np.argmax(objective)]

            error = np.max(np.abs(V_new - V_old)) #sup norm
            V_old = np.copy(V_new)

            count = count + 1

            if count > 1:
                ax.plot(y, V_old, color=plt.cm.jet(count/89), lw=1, alpha=0.6)

        ax.set_ylabel('Value', fontsize=12)
        ax.set_xlabel('Capital $k$', fontsize=12)
        ax.set_title('Value function iterations')
        plt.show()

        return policyFunction, V_new, count
