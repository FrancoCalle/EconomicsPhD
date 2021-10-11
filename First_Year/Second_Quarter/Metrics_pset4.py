import numpy as np
import matplotlib.pyplot as plt

R = np.array([1, 1, -1]).reshape(1,3)
I = np.eye(3,3)
A = np.array([[0.704, 0,0],[0, 1, -0.154],[0, 1, 0]])
B = np.array([[.144, 0], [0 , 0.206], [0,0]])

D = np.dot(R,A)
F = np.dot(R,B)


#Permanent Shock
ϕₚ = F + np.dot(np.dot(D, np.linalg.inv(I - A)), B)
ϕₜ = np.array([np.dot(np.dot(np.dot(D,
                        np.linalg.inv(I - A)),
                        np.linalg.matrix_power(A,i-1)), B).reshape(2,) for i in range(1,101)])


ϕ = ϕₚ-ϕₜ

# Impulse Response
first_shock, = plt.plot(ϕ[:,0], label = 'W1')
second_shock,  = plt.plot(ϕ[:,1], label = 'W2')
plt.legend(handles=[first_shock, second_shock])


## Kalman Filter Iteration:
#Initializer:
Σ0 = np.array([[0.2, 0.5, 0.8], [0.234, 0.532, 0.56],[0.1, 0.35, 0.66]])
eps = 0.0000001
diff = 10
ii = 0

Σₜ = np.copy(Σ0)

while (ii < 200) & (diff > eps):

        W = np.dot(np.dot(A,Σₜ), np.transpose(D)) + np.dot(B, np.transpose(F))
        Ω = np.dot(np.dot(D,Σₜ), np.transpose(D)) + np.dot(F, np.transpose(F))[0][0]
        XX = np.dot(np.dot(A, Σₜ), np.transpose(A)) + np.dot(B, np.transpose(B))

        Σ0 = np.copy(Σₜ)
        Σₜ = np.copy(XX - np.dot(W, np.transpose(W))/Ω)

        diff = np.abs(np.sum(Σ0 - Σₜ))
        ii +=1

# While Kalman Filter is:

W = np.dot(np.dot(A,Σₜ), np.transpose(D)) + np.dot(B, np.transpose(F))
κ = W/Ω


# Parameters for innovations representation

F_hat = np.sqrt(Ω)[0][0]
C = κ*F_hat


## Our contraction now is Σₜ

## Get impulse response of this new model:
#Permanent Shock
ϕₚ = F_hat + np.dot(np.dot(D, np.linalg.inv(I - A)), C)
ϕₜ = np.array([np.dot(np.dot(np.dot(D,
                        np.linalg.inv(I - A)),
                        np.linalg.matrix_power(A,i-1)), C).reshape(1,) for i in range(1,101)])

ϕ = ϕₚ-ϕₜ

# Impulse Response
plt.plot(ϕ[:,0])
