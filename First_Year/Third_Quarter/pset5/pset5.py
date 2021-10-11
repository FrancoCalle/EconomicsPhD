import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Static Parameters:
def setStaticParameters():

    nObservations = 1000
    nStates = 50
    interventionYearList = [1986, 1992, 1998, 2004]
    nGroups = len(interventionYearList)
    yearList = list(range(1980, 2011))
    delta = 0

    return nObservations, nStates, interventionYearList, nGroups, yearList, delta

# Random Allocation
def randomAllocation():

    #Random assignment to state:
    stateAllocation = list(np.random.randint(nStates, size=nObservations))

    #Random assignment of each state to treatment group:
    groupAllocation_i = np.random.randint(nGroups, size=nStates)

    #Look up state and assign treatment:
    groupAllocation = list(np.array(interventionYearList)[list(np.array(groupAllocation_i)[stateAllocation])])

    return stateAllocation, groupAllocation

# Set Parameters:
def setParameters(yearList,stateAllocation):

    #Time Fixed Effects:
    epsilonTimeFixedEffect = np.random.normal(0,1, len(yearList))
    #Individual fixed Effect:
    alpha_i_list = [np.random.normal(ss/5,1) for ss in stateAllocation]

    return epsilonTimeFixedEffect, alpha_i_list


# Generate Fake Data
def generateFakeData(nObservations, nStates, interventionYearList, nGroups, yearList, delta, groupAllocation, stateAllocation, epsilonTimeFixedEffect, alpha_i_list, seed = 84949):

    #Treatment effect:
    bbeta = 1

    #Immediate treamtment effect:
    np.random.seed(seed)

    data = {}


    for ii in range(nObservations):

        # Pick individual data:
        gg = groupAllocation[ii]     # Treatment year group

        #Define dictionary where we will store data:
        data[ii] = {}

        outcome_i = []
        tau_it_list = []
        epsilon_it_list = []
        alpha_t_list  = []

        for tt in yearList:

            # Obtain fixed effects parameters:
            alpha_i = alpha_i_list[ii]          # Individual fixed effect
            epsilon_it = np.random.normal(0, (.5)**2)   # No autocorrelation

            # Compute trend specific parameter:
            alpha_t = 0.1 * (tt - gg) + epsilonTimeFixedEffect[yearList.index(tt)]
            tau_it = bbeta*(tt - gg + 1) * (tt >= gg)

            # Generate Outcome:
            y_it = (2010 - gg) + alpha_i + alpha_t + tau_it + epsilon_it

            # Append outcomes:
            outcome_i.append(y_it)
            tau_it_list.append(tau_it)
            alpha_t_list.append(alpha_t)
            epsilon_it_list.append(epsilon_it)

        data[ii]['y_it'] = outcome_i
        data[ii]['tau_it'] = tau_it_list         #tau for each i and t
        data[ii]['alpha_t'] = alpha_t_list          #Outcome for each i and t
        data[ii]['epsilon_it'] = epsilon_it_list #Idiosincratic shock for each i and t

    return data


# set Dummies for TWFE regression:
def setDataForRegression(data):

    outcome = []
    year = []
    group = []
    id = []

    for ii in data:
        outcome += data[ii]['y_it']
        year += yearList
        group += [groupAllocation[ii] for xx in yearList]
        id += [ii for xx in yearList]

    DataFrame = pd.DataFrame({'id':id, 'outcome': outcome, 'year':year , 'group': group})
    DataFrame['event_year'] =  DataFrame['year']-DataFrame['group']
    DataFrame['d_event_year_l5'] =  DataFrame['event_year']<-5
    DataFrame['d_event_year_g5'] =  DataFrame['event_year']> 5

    # Dummy for event time between -5 and 5 inclusive:
    for dd in range(-5,6):
        if dd!=-1: DataFrame['event_year' + str(dd)] = (DataFrame['event_year'] == dd)

    # Dummy for each year:
    for yy in set(yearList):
        if yy != 2010: DataFrame['year_' + str(yy)] = (DataFrame['year'] == yy)

    # Dummy for each indivudual (individual Fixed Effect):
    for ii in set(id): DataFrame['id_' + str(ii)] = (DataFrame['id'] == ii)

    # Array with whole dataset
    XX = np.array(DataFrame.iloc[:,5:]).astype(int)
    Y = np.array(DataFrame.iloc[:,1])

    # Save column names:
    colnames  = list(DataFrame.iloc[:,5:].columns)

    return XX, Y, DataFrame, colnames


#OLS regression:
def olsRegression(XX, Y):

    bbeta = np.dot(np.linalg.inv(np.dot(XX.T,XX)), np.dot(XX.T,Y))

    return bbeta


#Monte Carlo simulations:
def monteCarloSimulations(nObservations,
                        nStates,
                        interventionYearList,
                        nGroups,
                        yearList,
                        delta,
                        groupAllocation,
                        stateAllocation,
                        epsilonTimeFixedEffect,
                        alpha_i_list):

    np.random.seed(0)
    seedList = np.random.randint(1,10**5,100)

    beta_s_list = []

    for seed in seedList:
        data_s = generateFakeData(nObservations, nStates, interventionYearList,
                                    nGroups, yearList, delta,
                                    groupAllocation, stateAllocation, epsilonTimeFixedEffect,
                                    alpha_i_list, seed)

        #Set data for regressions:
        XX, Y, _ = setDataForRegression(data_s)

        #Compute OLS regression:
        beta_s = olsRegression(XX, Y)


        beta_s_list.append(beta_s)

    return beta_s_list


#Estimate ATT(g,t) as in Callaway Santanna (2020):

def attCallawaySantanna(data):

    estimates = {}

    for g in range(3):#interventionYearList[:3]:

        gg = interventionYearList[g]
        estimates[gg] = {}

        for tt in yearList[1:24]: #iterate 5 periods over
            # First Difference:
            y_t_mean = DataFrame.loc[(DataFrame['year']==tt) &
                                        (DataFrame['group']==gg),'outcome'].mean() #Outcome when year t equals year g
            y_lag_mean = DataFrame.loc[(DataFrame['year']==gg-1) &
                                        (DataFrame['group']==gg),'outcome'].mean() #Outcome when year t equals year g

            d1 = y_t_mean - y_lag_mean

            # Second Difference:
            y_t_mean = DataFrame.loc[(DataFrame['year']==tt) &
                                        (DataFrame['group']==interventionYearList[-1]) &
                                        (DataFrame['event_year']<0),'outcome'].mean() #Outcome when year t equals year g

            y_lag_mean = DataFrame.loc[(DataFrame['year']==gg-1) &
                                        (DataFrame['group']==interventionYearList[-1]) &
                                        (DataFrame['event_year']<0),'outcome'].mean() #Outcome when year t equals year g

            d2 = y_t_mean - y_lag_mean

            # Compute difference:
            estimates[gg][tt] = d1-d2

    # Compile all estimates in Dataframe Format:
    dataframe = pd.DataFrame()

    for gg in estimates:
        coefficients = pd.DataFrame.from_dict(estimates[gg], orient='index').reset_index()
        coefficients['Group'] = gg
        dataframe = dataframe.append(coefficients)

    dataframe.rename(columns={"index": "Year", 0: "ATT"}, inplace= True)
    dataframe['event_year'] = dataframe['Year'] - dataframe['Group']

    return estimates, dataframe


# Aggregate parameters:
def aggregateCSestimates(dataframe):

    # Compute simple averages:
    ATT86 = dataframe.loc[(dataframe['Group'] == 1986) & (dataframe['event_year'] >= 0),'ATT'].mean()
    ATT92 = dataframe.loc[(dataframe['Group'] == 1992) & (dataframe['event_year'] >= 0),'ATT'].mean()
    ATT98 = dataframe.loc[(dataframe['Group'] == 1998) & (dataframe['event_year'] >= 0),'ATT'].mean()

    return ATT86, ATT92, ATT98


#-------------------------------------------------------------------------------
# Execute Stuff:
#-------------------------------------------------------------------------------

nObservations, nStates, interventionYearList, nGroups, yearList, delta =  setStaticParameters()

stateAllocation, groupAllocation = randomAllocation()

epsilonTimeFixedEffect, alpha_i_list = setParameters(yearList,stateAllocation)

data = generateFakeData(nObservations,
                        nStates,
                        interventionYearList,
                        nGroups,
                        yearList,
                        delta,
                        groupAllocation,
                        stateAllocation,
                        epsilonTimeFixedEffect,
                        alpha_i_list)

index_group1 = [i for i, x in enumerate(groupAllocation) if x == interventionYearList[0]]
index_group2 = [i for i, x in enumerate(groupAllocation) if x == interventionYearList[1]]
index_group3 = [i for i, x in enumerate(groupAllocation) if x == interventionYearList[2]]
index_group4 = [i for i, x in enumerate(groupAllocation) if x == interventionYearList[3]]

ATT_1 = np.array([data[ii]['tau_it'] for ii in index_group1])
ATT_2 = np.array([data[ii]['tau_it'] for ii in index_group2])
ATT_3 = np.array([data[ii]['tau_it'] for ii in index_group3])
ATT_4 = np.array([data[ii]['tau_it'] for ii in index_group4])

allyears = np.array(yearList)

# Part a: Plot ATT for each group in calendar time:
plt.plot(allyears, ATT_1.mean(0), label = 'Group: 1986')
plt.plot(allyears, ATT_2.mean(0), label = 'Group: 1992')
plt.plot(allyears, ATT_3.mean(0), label = 'Group: 1998')
plt.plot(allyears, ATT_4.mean(0), label = 'Group: 2004')
plt.legend()
plt.savefig('pset5/figure_a.png')
plt.show()

#Part B: plot ATT in "event time": Treatment effect is homogeneous.

plt.plot(allyears - interventionYearList[0], ATT_1.mean(0), label = 'Group: 1986')
plt.plot(allyears - interventionYearList[1], ATT_2.mean(0), label = 'Group: 1992')
plt.plot(allyears - interventionYearList[2], ATT_3.mean(0), label = 'Group: 1998')
plt.plot(allyears - interventionYearList[3], ATT_4.mean(0), label = 'Group: 2004')
plt.legend()
plt.savefig('pset5/figure_b.png')
plt.show()


# Part C: Yes I would ... blah blah
#----------------------------------


# Part D:
#--------

#First test that OLS is working as expected:
XX, Y, DataFrame, colnames = setDataForRegression(data)

#Compute OLS regression: (It's working well)
beta_s = olsRegression(XX, Y)

plt.hist(np.dot(XX,beta_s), bins =20)
plt.hist(Y, bins =20)

#Now do monte carlo simulations:
beta_s_list = monteCarloSimulations(nObservations,
                        nStates,
                        interventionYearList,
                        nGroups,
                        yearList,
                        delta,
                        groupAllocation,
                        stateAllocation,
                        epsilonTimeFixedEffect,
                        alpha_i_list)

# State column names of \gammas and index its location from colnames:
gammaNameList = ['event_year' + str(dd) for dd in range(-5,6) if dd != -1]
gammaIndex = [colnames.index(name) for name in gammaNameList]

beta_s_average = np.mean(np.array(beta_s_list),0)

gammas = beta_s_average[gammaIndex]
gamma_hat = np.insert(gammas,4,[0])

# Plot gamma estimates: (They are off indeed)
plt.scatter(range(-5,6), gamma_hat, label = 'γ estimates', color = 'blue')
plt.scatter(allyears - interventionYearList[1], ATT_2.mean(0), label = 'True treatment effect', color = 'orange') # Since the treatment effect is always the same we use group 2:
plt.ylim(-.5,7)
plt.xlim(-5.5,5.5)
plt.legend()
plt.savefig('pset5/figure_d.png')
plt.show()


# Part E: Is γ a measure of pre-trend?:
#------------------------------------

'''
Yes I would find pre trends since the slope is greater than zero meaning that the
slope for years prior to the policy in each group were positive. But in reality
we don't have pretrends in our specification
'''

# Part F: Discard TWFE and use Callaway & Sant'Anna (2020).
#----------------------------------------------------------

# Part G: Estimate ATT(g,t) based on f) and Plot estimates:
#----------------------------------------------------------

estimates, attCaSaDataFrame = attCallawaySantanna(data)

plt.scatter(attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1986,'Year'],attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1986,'ATT'], label = 'Group: 1986',color ="blue")
plt.scatter(attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1992,'Year'],attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1992,'ATT'], label = 'Group: 1992',color ="orange")
plt.scatter(attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1998,'Year'],attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1998,'ATT'], label = 'Group: 1998',color ="green")
plt.legend()
plt.savefig('pset5/figure_g.png')
plt.show() # SEXY!!!


# Part H: Obtain weighted average of ATT(g,t):
#---------------------------------------------

ATT86, ATT92, ATT98 = aggregateCSestimates(attCaSaDataFrame)

plt.scatter(attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1986,'Year'],attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1986,'ATT'], label = 'Group: 1986',color ="blue")
plt.scatter(attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1992,'Year'],attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1992,'ATT'], label = 'Group: 1992',color ="orange")
plt.scatter(attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1998,'Year'],attCaSaDataFrame.loc[attCaSaDataFrame['Group'] == 1998,'ATT'], label = 'Group: 1998',color ="green")
plt.axhline(y=ATT86, label = 'Avg Group: 1986',color ="blue", linestyle ="--")
plt.axhline(y=ATT92, label = 'Avg Group: 1992',color ="orange", linestyle ="--")
plt.axhline(y=ATT98, label = 'Avg Group: 1998',color ="green", linestyle ="--")
plt.legend()
plt.savefig('pset5/figure_h.png')
plt.show()
