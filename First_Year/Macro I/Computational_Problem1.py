#Import relevant packages
import numpy as np #gives you useful array operations
import matplotlib.pyplot as plt #this is for plotting
# %matplotlib inline


# Define parameters:

# bbeta = 0.99
# z = 1
# aalpha = 0.3
# ddelta = .1

parameters={
'beta': 0.99,
'z': 1,
'alpha': 0.3,
'delta':.1
}


gridSize = 500

# Define Grid:
y = np.linspace(0.01,15,gridSize) #state space
V_0 = [0 for k in y] #Initial guess forall states
V_old = np.copy(V_0)
V_new = np.copy(V_0)
policyFunction = np.zeros(V_old.shape)

################################################################################

#a) Value Function iteration:


def BellmanEquation(k, y, V0, parameters=parameters):

    bbeta = parameters['beta']
    z = parameters['z']
    aalpha = parameters['alpha']
    ddelta = parameters['delta']

    V = np.log(z*k**aalpha + (1-ddelta)*k - y) + bbeta*V0

    return V

def smartInitialGuessV0(y, parameters=parameters):

    bbeta = parameters['beta']
    z = parameters['z']
    aalpha = parameters['alpha']
    ddelta = parameters['delta']

    V_old = np.log(z*y**aalpha - ddelta*y)/(1-bbeta)

    return V_old

def steadyStateVariables(parameters):

    bbeta = parameters['beta']
    z = parameters['z']
    aalpha = parameters['alpha']
    ddelta = parameters['delta']

    #Actual Steady State of the model:
    steadyStateCapital = ((bbeta * z * aalpha)/(1- bbeta + ddelta*bbeta))**(1/(1-aalpha))
    steadyStateConsumption = z * steadyStateCapital**aalpha  - ddelta*steadyStateCapital

    return steadyStateCapital, steadyStateConsumption


def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return a.flat[idx]


def valueFunctionIteration(V_old, parameters=parameters):

    bbeta = parameters['beta']
    z = parameters['z']
    aalpha = parameters['alpha']
    ddelta = parameters['delta']

    fig, ax = plt.subplots()
    ax.plot(y, V_old, color=plt.cm.jet(0), lw=1, alpha=0.6)

    error = 9999
    tol = 1e-5
    count = 0

    while error>tol:

        V_new = np.ones(gridSize)

        for i in range(gridSize):
            k = y[i]
            upperBoundY = (z*k**aalpha + (1-ddelta)*k)
            objective = [BellmanEquation(k,y[j],V_old[j]) if y[j] <= upperBoundY else -99999 for j in range(gridSize)]
            V_new[i] = np.max(objective)
            policyFunction[i] = y[np.argmax(objective)]

        error = np.max(np.abs(V_new - V_old)) #sup norm
        V_old = np.copy(V_new)

        count = count + 1

        ax.plot(y, V_old, color=plt.cm.jet(count/89), lw=1, alpha=0.6)

    ax.set_ylabel('Value', fontsize=12)
    ax.set_xlabel('Capital $k$', fontsize=12)
    ax.set_title('Value function iterations')
    plt.show()

    return policyFunction, V_new, count

parameters['beta'] = 0.95
policyFunction, V_new, count = valueFunctionIteration(V_old)

print('Number of iterations:', count)


# Plot Value Function, All iterations


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
c_ss_numerical = parameters['z']*k**parameters['alpha'] + (1-parameters['delta'])*k_ss_numerical - k_ss_numerical
print(" Capital Steady State", round(k_ss_numerical,3),"\n Consumption Steady State:", round(c_ss_numerical,3))




#Actual Steady State of the model:
steadyStateCapital, steadyStateConsumption = steadyStateVariables(parameters)
print(" Capital Steady State", round(steadyStateCapital,3),"\n Consumption Steady State:", round(steadyStateConsumption,3))


parameters['beta'] = 0.1
policyFunction, V_new, count = valueFunctionIteration(V_old)
print(count)

# Plot latest value function:
fig, ax2 = plt.subplots()
ax2.plot(y , V_new, color=plt.cm.jet(count/500), linewidth=2.5)
ax2.set_ylabel('value', fontsize=12)
ax2.set_xlabel('Capital $k$', fontsize=12)
ax2.set_title('Value function iterations')


# Now let's use a new initial values for v0:

parameters['beta'] = 0.99
V_old = smartInitialGuessV0(y)
policyFunction, V_new, count = valueFunctionIteration(V_old)

# Question 3:
parameters['z'] = 1
V_old = smartInitialGuessV0(y)
policyFunctionZ1, V_new, count = valueFunctionIteration(V_old)

parameters['z'] = 2
V_old = smartInitialGuessV0(y)
policyFunctionZ2, V_new, count = valueFunctionIteration(V_old)


parameters['z'] = 1
steadyStateCapital, steadyStateConsumption = steadyStateVariables(parameters)

k = find_nearest(y,steadyStateCapital)
error = 10
k_transition = [steadyStateCapital]*4

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
