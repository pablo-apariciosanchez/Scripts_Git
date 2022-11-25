# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 09:03:36 2020

"""

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer
import time


df =  pd.read_csv('C:/Users/Usuario/Desktop/WOTAN_parallel_3.0/DatasetB-calculated_with_y.csv', sep=';')
df.head()

ns = time.time()


columns = df[['SMILES', 'y']]

df_features = df.iloc[:, 2:]
df_features.head()
df.dtypes
df_features.dtypes

column_names = list(df_features.columns)

df_nan = df_features.replace(to_replace = r'[a-zA-Z]', value = np.nan, regex = True).replace(np.inf, np.nan)

# print(df_nan)

to_cast = list(df_nan.select_dtypes(include=[object]).columns)
# print(to_cast)

df_nan[to_cast] = df_nan[to_cast].astype(dtype=np.float64)
# df_nan.dtypes

# df_dropped = df_nan.drop(to_cast)

df_nan[df_nan <= -1e38] = np.nan
df_nan[df_nan >= 1e38] = np.nan

imputer = KNNImputer(missing_values = np.nan, n_neighbors=3, weights="uniform")
result = imputer.fit_transform(df_nan)
# print(result)

result_df = pd.DataFrame(result, columns = column_names)
# print(result_df)

ns_final = time.time()

print('Time: ', ns_final-ns, 'seconds')


df_nan2 = columns.join(result_df)
df_nan2.head()


df_nan2.to_csv('C:/Users/Usuario/Desktop/Scripts/data/AMES_calculated_with_y_knn.csv', sep=';', index = False)


###############################################################################
###############################################################################

dfa = pd.read_csv('C:/Users/Usuario/Desktop/DILI_model/DILI_all_sanitized_cleaned_for_processandmodel_calculated_with_y_new_mordred_knn.csv', sep=';')
dfa.head()

dfa_feat = dfa.drop(columns = ['SMILES', 'y'], axis = 1)

dfb = pd.read_csv('C:/Users/Usuario/Desktop/DILI_model/DILI_all_sanitized_cleaned_for_processandmodel_calculated_with_y_wotan_knn.csv', sep=';')

df_final = dfb.join(dfa_feat)
df_final.head()

df_final.to_csv('C:/Users/Usuario/Desktop/DILI_model/DILI_all_sanitized_cleaned_for_processandmodel_calculated_with_y_ALL_knn.csv', sep=';', index = False)


