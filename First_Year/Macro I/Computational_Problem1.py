#Import relevant packages
import numpy as np #gives you useful array operations
import matplotlib.pyplot as plt #this is for plotting
from value_function_iteration import dynamicProgramingModule
# %matplotlib inline


# Define Grid:
y = np.linspace(0.01,15,500) #state space
V_0 = [0 for k in y] #Initial guess forall states
V_old = np.copy(V_0)
V_new = np.copy(V_0)
policyFunction = np.zeros(V_old.shape)

################################################################################


#a) Value Function iteration:

bbeta = 0.99
z = 1
dpm = dynamicProgramingModule(bbeta, z)
policyFunction, V_new, count = dpm.valueFunctionIteration(V_old,y)

print('Number of iterations:', count)


# Plot Value Function, Final iteration
fig, ax2 = plt.subplots()
ax2.plot(y , V_new, color=plt.cm.jet(count/500), linewidth=2.5)
ax2.set_ylabel('value', fontsize=12)
ax2.set_xlabel('Capital $k$', fontsize=12)
ax2.set_title('Value function iterations')


# Plot Policy Rule, Final iteration
fig, ax3 = plt.subplots()
ax3.plot(y , policyFunction, color=plt.cm.jet(count/500), linewidth=2.5)
ax3.set_ylabel('Policy Function $y$', fontsize=12)
ax3.set_xlabel('Capital $k$', fontsize=12)
ax3.set_title('Value function iterations')
plt.show()


#Steady State Numerical Approximation:

# Start with k = k[1], then apply g(k) to get k' and then continue iteratively:
k = y[0]
error = 10
k_transition = []

count = 1
while (error > 1e-5) & (count < 200):
    k_new = policyFunction[k==y][0]
    error = np.abs(k_new-k)
    k_transition.append(k_new)
    k = k_new
    count += 1


fig, ax4 = plt.subplots()
ax4.plot(k_transition, color=plt.cm.jet(count/500), linewidth=2.5)
ax4.set_ylabel('Policy Function $y$', fontsize=12)
ax4.set_xlabel('Time', fontsize=12)
ax4.set_title('Steady State')
plt.show()

#Steady State approximation:
k_ss_numerical = k_new
c_ss_numerical = dpm.z*k**dpm.aalpha + (1-dpm.ddelta)*k_ss_numerical - k_ss_numerical
print(" Capital Steady State", round(k_ss_numerical,3),"\n Consumption Steady State:", round(c_ss_numerical,3))


#Actual Steady State of the model:
steadyStateCapital, steadyStateConsumption = dpm.steadyStateVariables()
print(" Capital Steady State", round(steadyStateCapital,3),"\n Consumption Steady State:", round(steadyStateConsumption,3))


#Evaluate with beta = 0.1
dpm1 = dynamicProgramingModule(0.1, 1)
policyFunction, V_new, count = dpm1.valueFunctionIteration(V_old,y)
print(count)

# Plot latest value function:
fig, ax2 = plt.subplots()
ax2.plot(y , V_new, color=plt.cm.jet(count/500), linewidth=2.5)
ax2.set_ylabel('value', fontsize=12)
ax2.set_xlabel('Capital $k$', fontsize=12)
ax2.set_title('Value function iterations')


# Now let's use a new initial values for v0:
V_old = dpm.smartInitialGuessV0(y)
policyFunction, V_new, count = dpm.valueFunctionIteration(V_old,y)

# Question 3:
dpmZ1 = dynamicProgramingModule(0.99, 1)
V_old = dpmZ1.smartInitialGuessV0(y)
policyFunctionZ1, V_new, count = dpmZ1.valueFunctionIteration(V_old,y)
steadyStateCapitalZ1, steadyStateConsumptionZ1 = dpmZ1.steadyStateVariables()

dpm3 = dynamicProgramingModule(0.99, 2)
V_old = dpm3.smartInitialGuessV0(y)
policyFunctionZ2, V_new, count = dpm3.valueFunctionIteration(V_old,y)


# Plot productivity shock simulation
k = dpm3.find_nearest(y,steadyStateCapitalZ1)
error = 10
k_transition = [steadyStateCapitalZ1]*4

for ii in range(29):
    k_new = policyFunctionZ2[k==y][0]
    error = np.abs(k_new-k)
    k_transition.append(k_new)
    k = k_new
    count += 1


fig, ax4 = plt.subplots()
ax4.plot(k_transition, color=plt.cm.jet(count/500), linewidth=2.5)
ax4.set_ylabel('Capital $k$', fontsize=12)
ax4.set_xlabel('Time', fontsize=12)
ax4.set_title('Transition to Steady State')
plt.show()
