# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 15:00:13 2021

"""


import pyreadstat
import pandas as pd
import pickle
import os
os.chdir('C:/Users/proto/Desktop/Pablo/ProtoADME/ProtoADME_v1_0/properties/endpoints/Half_life')



model = pickle.load(open('Half_life.sav', 'rb'))

model.get_params()
print(model.best_params_) 
print(model.scoring) 


