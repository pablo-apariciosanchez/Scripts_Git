# -*- coding: utf-8 -*-
"""
Created on Thu Jun 9 09:14:24 2022

"""


import os
# os.chdir('F:/ProtoQSAR/20200318/some_miscelaneous_scripts/AD')
os.chdir('C:/Users/Usuario/Desktop/Scripts/Data/')

import pandas as pd
import numpy as np

from sklearn.preprocessing import StandardScaler, MinMaxScaler
import matplotlib.pyplot as plt
from rdkit.Chem import AllChem
from rdkit import Chem

###############################################################################

# train_data = pd.read_csv('melting_point-descriptors-train.txt', sep='\t', header=0, encoding='utf-8')
train_data = pd.read_csv('dataset-descriptors-train.txt', sep='\t', header=0, encoding='utf-8')

train_data.head()


###############################################################################
################################# LEVERAGES ###################################
###############################################################################

# select descriptors part of df
descriptors_train = train_data.iloc[:,2:]
descriptors_train = train_data.iloc[:, 1:-2]
descriptors_train.head()

# scale 

scaler = StandardScaler()

scaler.fit(descriptors_train)

X = scaler.transform(descriptors_train)

# X = np.array(descriptors_train)

# h value for each descriptor row in dataset is defined as hii = xi T(XTX)–1 xi

## calculate (XTX)–1

tranpX = np.transpose(X)

xtx = np.dot(tranpX,X)

invxtx = np.linalg.inv(xtx)

## calculate h value for each row:

leverages = []

for row in X:
    transposed_row = np.transpose(row)
    leverage = np.dot(np.dot(row, invxtx), transposed_row)
    leverages.append(leverage)


# calculate leverage warning:
p = X.shape[1]
n = X.shape[0]
warning_leverage = 3*(p/n)



# plot
warning_leverage_line = [warning_leverage] * n
x_ax = range(len(list(train_data['SMILES'])))


outliers_train = []


for i,lev in enumerate(leverages):
    if lev > warning_leverage:
        outliers_train.append(i)
        print(i,lev)


for i,row in train_data.iterrows():
    if i in outliers_train:
        print(i, row[0])




plt.scatter(x_ax, leverages, s=5, color="blue", label="h")
plt.plot(x_ax, warning_leverage_line, lw=1.5, color="red", label="h*")
plt.xlabel('sample', fontsize=10)
plt.ylabel('leverage', fontsize=10)

# for x,y in zip(x_ax,leverages):

#     label = x

#     plt.annotate(label, # this is the text
#                   (x,y), # this is the point to label
#                   textcoords="offset points", # how to position the text
#                   xytext=(0,10),
#                   fontsize=7,# distance from text to points (x,y)
#                   ha='center') 


# plt.legend()
plt.show()


#%%

#################### Standardized residuals vs. leverages (William's plot) ###############
"""
@author: Pablo Aparicio
"""

# train_data = pd.read_csv('C:/Users/Usuario/Desktop/Caco2_test/Test5/Caco2-descriptors-train.txt', sep='\t', header=0, encoding='utf-8')
# train_data.head()

import math

y_obs = train_data['y']
y_pred = train_data['predicted']

# len(y_obs)
# len(y_pred)
    
residue = y_obs - y_pred

residues = residue.tolist()

# media = sum(residues)/len(y_obs)
media = np.mean(residues)

# var = sum((residue-media)**2)/n
var = np.var(residues)

sqrt_var = math.sqrt(var)

se_regs =[]

for lev in leverages:
   se_reg = sqrt_var*((1-lev)**0.5)
   se_regs.append(se_reg)

studentized_residuals = residue/se_regs

st_res = studentized_residuals.to_list()


warning_leverage = 3*(p/n)

# plot
warning_leverage_line = [warning_leverage] *n


residuos_down= [-3] * n
residuos_up= [3] * n

plt.scatter(leverages, st_res, s=5, color="blue", label="h")

plt.plot(leverages, residuos_down, lw=1.5, color="red")
plt.plot(leverages, residuos_up, lw=1.5, color="red")
# plt.plot(leverages, residuos_down, lw=1.5, color="red")
plt.plot(warning_leverage_line, st_res, lw=1.5, color="green")


plt.xlabel('Leverage', fontsize=10)
plt.ylabel('Standarized Residues', fontsize=10) 

plt.show()


##############################################################################
#%%

# test_data = pd.read_csv('files/{}-descriptors-test.txt'.format(model), sep='\t', header=0, encoding='utf-8')
test_data = pd.read_csv('C:/Users/Usuario/Desktop/DILI_model/DILI-descriptors-test.txt', sep='\t', header=0, encoding='utf-8')

descriptors_test = test_data.iloc[:, 1:-2]

print(descriptors_test.dtypes)
Y = scaler.transform(descriptors_test)


leverages_test = []

for row in Y:
    transposed_row = np.transpose(row)
    leverage = np.dot(np.dot(row, invxtx), transposed_row)
    leverages_test.append(leverage)


outliers_test = []

for i,lev in enumerate(leverages_test):
    if lev > warning_leverage:
        outliers_test.append(i)
        print(i,lev)


for i,row in test_data.iterrows():
    if i in outliers_test:
        print(i, row[0])



n_test = Y.shape[0]

warning_leverage_line2 = [warning_leverage] * n_test

x_ax_test = range(len(list(test_data['SMILES'])))


plt.scatter(x_ax_test, leverages_test, s=5, color="red", label="h")
plt.plot(x_ax_test, warning_leverage_line2, lw=1.5, color="blue", label="h*")
plt.xlabel('sample', fontsize=10)
plt.ylabel('leverage', fontsize=10)

# for x,y in zip(x_ax,leverages):

#     label = x

#     plt.annotate(label, # this is the text
#                  (x,y), # this is the point to label
#                  textcoords="offset points", # how to position the text
#                  xytext=(0,10),
#                  fontsize=7,# distance from text to points (x,y)
#                  ha='center')


# plt.legend()
plt.show()

