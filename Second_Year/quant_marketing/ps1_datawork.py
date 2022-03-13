import pandas as pd 
import numpy as np
import copy
import warnings
warnings.filterwarnings('ignore')


df = pd.read_csv("C:/Users/franc/Dropbox/Franco Econ Phd/2 Second Year/Winter/Quantitative Marketing/ps1/pbout_2.csv", delimiter=',', index_col=False)

#%%

price_data = (df.groupby(["date", "choice"])['price']
                .mean()
                .reset_index()
                .sort_values("date")
                )

dates = np.unique(price_data.date)

date_price_dictionary = {}

for ii in range(np.size(dates)):

    dd = dates[ii]

    date_price = (price_data.loc[price_data.date == dd,['choice', 'price']]
                                    .groupby('choice')
                                    .price
                                    .mean()
                                    )
    
    jj = copy.copy(ii)

    while (len(list(date_price.index)) < 9) & (jj < len(dates)-1):

        choice_list_in_date = list(date_price.index)
        
        jj = jj + 1

        dd_next = dates[jj]

        date_price = date_price.append(
                                    (price_data.loc[(price_data.date == dd_next) & ~(price_data.choice.isin(choice_list_in_date)),['choice', 'price']]
                                    .groupby('choice')
                                    .price
                                    .mean()
                                    )
                                    )
    
    if len(date_price.index) == 9:

        date_price_dictionary[dd] = date_price.sort_index().to_dict()

    else:

        date_price_dictionary[dd] = date_price_dictionary[dates[ii-1]]



len(list(date_price_dictionary))

price_data_df = pd.DataFrame(date_price_dictionary).T

price_data_df.columns = ["price_" + str(cc) for cc in price_data_df.columns]

price_data_df = price_data_df.reset_index().rename(columns={'index':'date'})

price_data_df.to_csv("C:/Users/franc/Dropbox/Franco Econ Phd/2 Second Year/Winter/Quantitative Marketing/ps1/price_by_date.csv", index=False)

#%%

df = df.merge(price_data_df, on="date").drop(columns='Unnamed: 0')

choice_cols = pd.get_dummies(df.choice)

choice_cols.columns = ["choice_" + str(cc) for cc in choice_cols.columns]

df = pd.concat([df, choice_cols], axis=1)

df.to_csv("C:/Users/franc/Dropbox/Franco Econ Phd/2 Second Year/Winter/Quantitative Marketing/ps1/pbout_final.csv", index=False)

