import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

df = pd.read_csv('ps3_auction.csv')
# df["optim_rn_Bid6"] = df.BidC6 * 5/6
# df["optim_rn_Bid3"] = df.BidC3 * 2/3


# Scatter plot (A)
plt.scatter(df.Value, df.BidC3, marker="+", color='gray')
plt.plot(df.Value, df.Value*(2/3), color='black')
plt.xlim(0,30)
plt.ylim(0,30)
plt.xlabel("Valuation")
plt.ylabel("Bid")


# Scatter plot (A)
plt.scatter(df.Value, df.BidC6, marker="+", color='gray')
plt.plot(df.Value, df.Value*(5/6), color='black')
plt.xlim(0,30)
plt.ylim(0,30)
plt.xlabel("Valuation")
plt.ylabel("Bid")



# Scatter plot (A)
df_ex_3 = df.loc[df.Experiment==3,:]
df_ex_4 = df.loc[df.Experiment==4,:]
df_ex_5 = df.loc[df.Experiment==5,:]

grid = np.array((range(0,3001,1)))/100

# Number of bidders: 3
empirical_cdf_3_nbid3 = np.array([np.mean(df_ex_3.BidC3 < x) for x in grid])
empirical_cdf_4_nbid3 = np.array([np.mean(df_ex_4.BidC3 < x) for x in grid])
empirical_cdf_5_nbid3 = np.array([np.mean(df_ex_5.BidC3 < x) for x in grid])

plt.plot(grid,empirical_cdf_3_nbid3, label="Experiment 3")
plt.plot(grid,empirical_cdf_4_nbid3, label="Experiment 4")
plt.plot(grid,empirical_cdf_5_nbid3, label="Experiment 5")
plt.xlim(0,25)
plt.ylim(0,1)
plt.xlabel("Bid")
plt.ylabel("F(Bid)")
plt.savefig("q2_cdf_nbidders_3.pdf")


# Number of bidders: 6
empirical_cdf_3_nbid6 = np.array([np.mean(df_ex_3.BidC6 < x) for x in grid])
empirical_cdf_4_nbid6 = np.array([np.mean(df_ex_4.BidC6 < x) for x in grid])
empirical_cdf_5_nbid6 = np.array([np.mean(df_ex_5.BidC6 < x) for x in grid])

plt.plot(grid,empirical_cdf_3_nbid6, label="Experiment 3")
plt.plot(grid,empirical_cdf_4_nbid6, label="Experiment 4")
plt.plot(grid,empirical_cdf_5_nbid6, label="Experiment 5")
plt.xlim(0,25)
plt.ylim(0,1)
plt.xlabel("Bid")
plt.ylabel("F(Bid)")
plt.savefig("q2_cdf_nbidders_6.pdf")



# Question 3: Empirical CDF to compute expected profit and optimization error:

# Expected profit: F(bid) * (v - bid)

# N Players: 6 

empirical_cdf_nbid3 = np.array([np.mean(df.BidC3 <= x) for x in grid])

empirical_cdf_nbid6 = np.array([np.mean(df.BidC6 <= x) for x in grid])


def obtain_optimization_error(df, empirical_cdf, N):

    expected_profit_list = []

    optimal_bid_list = []

    if N==6:
        bid_list = df.BidC6
    else:
        bid_list = df.BidC3


    for bb, vv in zip(bid_list, df.Value):
        
        F_bid = empirical_cdf[(np.abs(grid - bb)).argmin()]

        expected_profit = (vv - bb) * F_bid**(N-1)

        expected_profit_list.append(expected_profit)

        optimal_profit = np.max((vv - grid) * empirical_cdf**(N-1))                            
        
        # optimal_profit = grid[np.argmax((vv - grid) * empirical_cdf**(N-1))]                            

        optimal_bid_list.append(optimal_profit)

    optimal_bid_list = np.array(optimal_bid_list)

    expected_profit_list = np.array(expected_profit_list)
    
    optimization_error = optimal_bid_list - expected_profit_list

    return expected_profit_list, optimal_bid_list, optimization_error

# Q 3: 

# Number of bidders = 6:
exp_profit_list_nb6, opt_profit_list_nb6, opt_error_nb6 = obtain_optimization_error(df, empirical_cdf_nbid6, N = 6)

# Number of bidders = 3:
exp_profit_list_nb3, opt_profit_list_nb3, opt_error_nb3 = obtain_optimization_error(df, empirical_cdf_nbid3, N = 3)


# Q 4: Generate Table 1:

# Expected Profit under:

# N = 3
exp_profit_3 = [round(np.mean(exp_profit_list_nb3),2), round(np.percentile(exp_profit_list_nb3,25),2), 
                round(np.percentile(exp_profit_list_nb3,50),2), round(np.percentile(exp_profit_list_nb3,75),2)]

# N = 6
exp_profit_6 = [round(np.mean(exp_profit_list_nb6),2), round(np.percentile(exp_profit_list_nb6,25),2), 
                round(np.percentile(exp_profit_list_nb6,50),2), round(np.percentile(exp_profit_list_nb6,75),2)]


# Optimization Error under:

# N = 3
opt_error_3 = [round(np.mean(opt_error_nb3),2), round(np.percentile(opt_error_nb3,25),2), 
                round(np.percentile(opt_error_nb3,50),2), round(np.percentile(opt_error_nb3,75),2)]

# N = 6
opt_error_6 = [round(np.mean(opt_error_nb6),2), round(np.percentile(opt_error_nb6,25),2), 
                round(np.percentile(opt_error_nb6,50),2), round(np.percentile(opt_error_nb6,75),2)]




# Computations Section III:

def silverman_badwidth(x_list,df):
            
    n = df.shape[0]
    q75, q25 = np.percentile(x_list, [75 ,25])
    iqr = np.abs(q75 - q25)
    std_hat = np.std(x_list)
    h = .9 * min(std_hat, iqr/1.34) * n **(-.2)

    return h


def normal_kernel(x, df, N=3):
    
    if N == 3:
        x_list = df.BidC3
    else:
        x_list = df.BidC6

    h = silverman_badwidth(x_list, df)

    n = len(x_list)

    eps = np.array(x_list - x)/h

    den_i = (1/(n*h)) * 1/np.sqrt(2 * np.pi) * np.sum(np.exp(-1/2 * eps ** 2))

    return den_i 


bid_c3 = np.array(df.BidC3)
bid_c6 = np.array(df.BidC6)
empirical_pdf_n3 = np.array([normal_kernel(x, df, 3) for x in bid_c3])
empirical_pdf_n6 = np.array([normal_kernel(x, df, 6) for x in bid_c6])



plt.scatter(df.BidC3,empirical_pdf_n3, marker="+", color="gray")
plt.savefig("density_n3.pdf")

plt.scatter(df.BidC6,empirical_pdf_n6, marker="+", color="gray")
plt.savefig("density_n6.pdf")

bid_c3_sorted = bid_c3[np.argsort(bid_c3)]
epdf_n3_sorted = empirical_pdf_n3[np.argsort(bid_c3)]
ecdf_n3_sorted = np.cumsum(epdf_n3_sorted)
ecdf_n3_sorted = ecdf_n3_sorted/max(ecdf_n3_sorted)

bid_c6_sorted = bid_c6[np.argsort(bid_c6)]
epdf_n6_sorted = empirical_pdf_n6[np.argsort(bid_c6)]
ecdf_n6_sorted = np.cumsum(epdf_n6_sorted)
ecdf_n6_sorted = ecdf_n6_sorted/max(ecdf_n6_sorted)

plt.plot(bid_c3_sorted,ecdf_n3_sorted, color="gray")
plt.plot(bid_c6_sorted,ecdf_n6_sorted, color="gray")


# Now obtain valuations based on equation 9:

v_n6 = bid_c6_sorted + .20 * ecdf_n6_sorted/epdf_n6_sorted
plt.hist(v_n6, bins=60, color='gray', edgecolor="black")
plt.xlim(0,40)
plt.savefig("estimated_valuation_n6.pdf")


v_n3 = bid_c3_sorted + .50 * ecdf_n3_sorted/epdf_n3_sorted
plt.hist(v_n3, bins=70, color='gray', edgecolor="black")
plt.xlim(0,60)
plt.ylim(0,30)
plt.savefig("estimated_valuation_n3.pdf")


true_valuation_c3_sorted = np.array(df.Value)[np.argsort(bid_c3)]
true_valuation_c6_sorted = np.array(df.Value)[np.argsort(bid_c6)]

plt.scatter(true_valuation_c3_sorted, v_n3, marker="+", color="gray")
plt.plot(true_valuation_c3_sorted, true_valuation_c3_sorted, color="black")
plt.ylim(-2,80)
plt.savefig("true_vs_estimated_valuation_n3.pdf")

plt.scatter(true_valuation_c6_sorted, v_n6, marker="+", color="gray")
plt.plot(true_valuation_c6_sorted, true_valuation_c6_sorted, color="black")
plt.ylim(-2,60)
plt.savefig("true_vs_estimated_valuation_n6.pdf")

# L 1 Norm:
np.mean(abs(true_valuation_c6_sorted -  v_n6))
np.mean(abs(true_valuation_c3_sorted -  v_n3))

# L 2 Norm:
np.sqrt(1/len(v_n3)*np.sum((true_valuation_c3_sorted -  v_n3)**2))
np.sqrt(1/len(v_n6)*np.sum((true_valuation_c6_sorted -  v_n6)**2))


# SECTION IV:
# Risk Aversion: Part 1:
ptile = list(range(1,101))
ptile = list(range(5,96))
ptile = list(range(25,76))

def estimate_risk_aversion(ptile):

    dif_bb_list = []
    dif_v_hat_list = []

    for pp in ptile:
        
        bb_c3 = np.percentile(bid_c3,pp)
        pdf_bb = normal_kernel(bb_c3, df, 3)
        cdf_bb = np.mean(bid_c3<=bb_c3)
        v_hat_n3 = (1/2)*cdf_bb/pdf_bb

        bb_c6 = np.percentile(bid_c6,pp)
        pdf_bb = normal_kernel(bb_c6, df, 6)
        cdf_bb = np.mean(bid_c6<=bb_c6)
        v_hat_n6 = (1/5)*cdf_bb/pdf_bb

        dif_bb = bb_c3 - bb_c6 
        dif_v_hat = v_hat_n6 - v_hat_n3

        dif_bb_list.append(dif_bb)
        dif_v_hat_list.append(dif_v_hat)


    dif_bb_list = np.array(dif_bb_list)
    dif_v_hat_list = np.array(dif_v_hat_list)

    # OLS without intercept:
    # theta = sum(dif_v_hat_list**2)**(-1)*sum(dif_bb_list*dif_v_hat_list)
    est=sm.OLS(dif_bb_list, dif_v_hat_list)
    est = est.fit()
    theta = est.params[0]
    
    return theta, est

# TODO: Add 95 percent confidence intervals...

# Full Sample:
theta_full, est_full = estimate_risk_aversion(list(range(0,101)))
# OLS Trimming 5th tails:
theta_p5, est_p5 = estimate_risk_aversion(list(range(5,96)))
# OLS Trimming 25th tails:
theta_p25, est_p25 = estimate_risk_aversion(list(range(25,76)))

est_full.summary()
est_p5.summary()
est_p25.summary()

# Estimated Valuation distribution:
v_n3_risk_aversion = bid_c3_sorted + theta_p25 * .50 * ecdf_n3_sorted/epdf_n3_sorted
plt.hist(v_n3_risk_aversion, bins=40, color='gray', edgecolor="black")
plt.xlim(0,40)
plt.ylim(0,30)
plt.savefig("estimated_valuation_cra_n3.pdf")

v_n6_risk_aversion = bid_c6_sorted + theta_p25 * .20 * ecdf_n6_sorted/epdf_n6_sorted
plt.hist(v_n6_risk_aversion, bins=40, color='gray', edgecolor="black")
plt.xlim(0,35)
plt.ylim(0,25)
plt.savefig("estimated_valuation_cra_n6.pdf")


# Scatter plots:
plt.scatter(true_valuation_c3_sorted, v_n3_risk_aversion, marker="+", color="gray")
plt.plot(true_valuation_c3_sorted, true_valuation_c3_sorted, color="black")
plt.ylim(-2,32)
plt.savefig("true_vs_estimated_valuation_cra_n3.pdf")

plt.scatter(true_valuation_c6_sorted, v_n6_risk_aversion, marker="+", color="gray")
plt.plot(true_valuation_c6_sorted, true_valuation_c6_sorted, color="black")
plt.ylim(-2,32)
plt.savefig("true_vs_estimated_valuation_cra_n6.pdf")


# L 1 Norm:
np.mean(abs(true_valuation_c3_sorted -  v_n3_risk_aversion))
np.mean(abs(true_valuation_c6_sorted -  v_n6_risk_aversion))


# L 2 Norm:
np.sqrt(1/len(v_n3_risk_aversion)*np.sum((true_valuation_c3_sorted -  v_n3_risk_aversion)**2))
np.sqrt(1/len(v_n6_risk_aversion)*np.sum((true_valuation_c6_sorted -  v_n6_risk_aversion)**2))

