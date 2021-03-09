import numpy as np
import matplotlib.pyplot as plt

βz = np.exp(-.014)
σz = np.array([.011, 0.025])

βc = 1
σc = np.array([0.477, 0])

nper = 400

#Permanent Shock
ϕₚ = σc + βc * (1 - βz)**(-1) * σz
ϕₜ = np.array([βc * (1 - βz)**(-1) * βz**(i-1) * σz  for i in range(1,nper)])

ϕ = ϕₚ-ϕₜ

# Impulse Response
first_shock, = plt.plot(ϕ[:,0], label = 'W1')
second_shock,  = plt.plot(ϕ[:,1], label = 'W2')
plt.legend(handles=[first_shock, second_shock])
plt.savefig("IRF_Consumption_Model.png")


#Part b:

ϕₚ
ϕₒ = [1, -ϕₚ[0]/ϕₚ[1]]    # Ortogonal shock

#Ortogonalization matrix
Ω=np.array([ϕₚ/np.linalg.norm(ϕₚ), ϕₒ/np.linalg.norm(ϕₒ)])
φ = np.dot(ϕ,np.linalg.inv(Ω)) #Ortogonailzed shocks


first_shock, = plt.plot(φ[:,0], label = 'W1')
plt.legend(handles=[first_shock])
plt.savefig("Shock for [1 0]")

second_shock, = plt.plot(φ[:,1], label = 'W2')
plt.legend(handles=[second_shock])
plt.savefig("Shock for [0 1]")

first_shock, = plt.plot(φ[:400,0], label = 'W1')
second_shock,  = plt.plot(φ[:400,1], label = 'W2')
plt.legend(handles=[first_shock, second_shock])



# Risk Price:
γ_list = [1, 2, 4, 8]
λ = np.exp(-0.002)

σs = [np.dot(γ * σc - (1-γ)*λ*σz/(1-λ*βz),Ω) for γ in γ_list]
