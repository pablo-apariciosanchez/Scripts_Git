# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 15:26:56 2022

"""

import pandas as pd
import requests

df = pd.read_csv(r'C:\Users\proto\Desktop\Pablo\PAS\ADMET_properties\4.Excretion\Clearance\Renal CL\RenalCLdataset.csv', sep=',')


identifiers  = df.name


smiles_df = pd.DataFrame(columns = ['Name', 'Smiles'])

for x in identifiers :
    try:
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' + x + '/property/CanonicalSMILES/TXT'
#         remove new line character with rstrip
        smiles = requests.get(url).text.rstrip()
        if('NotFound' in smiles):
            print(x, " not found")
        else: 
            smiles_df = smiles_df.append({'Name' : x, 'Smiles' : smiles}, ignore_index = True)
    except: 
        print("boo ", x)
print(smiles_df)

smiles_df.to_csv(r'C:\Users\proto\Desktop\Pablo\PAS\ADMET_properties\4.Excretion\Clearance\Renal CL\RenalCLdataset_withSMILES.csv', sep=';', index=False)
