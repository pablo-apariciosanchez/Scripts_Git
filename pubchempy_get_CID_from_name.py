# -*- coding: utf-8 -*-
"""
Created on Fri May 27 13:23:37 2022

"""

import pubchempy as pcp
import os
import pandas as pd

# PATH = "F:/RITA/TESIS/datos_para_tesis/OBrien_HepG2_HCS"
# file = "OBrien_SMILES_amano_effectsIncrDecr.xlsx"
# col_name = "Drug_name"
# PATH = "F:\RITA\TESIS\datos_para_tesis\DILI_classification_databases\DILIst_binary_DILI_database_1200compounds"
# file = "DILIst_compounds.xlsx"
# col_name = "CompoundName"
# PATH = "F:\RITA\TESIS\datos_para_tesis\DILI_classification_databases\DILIrank"
# file = "DILIrank_compounds.xlsx"
# col_name = "Compound Name"
# PATH = "F:/RITA/TESIS/datos_para_tesis/LA_FE/Laia_combine_w_Tolosa2012"
# file = "TolosaNew_compounds_bin5props_3h24h.xlsx"
# col_name = "Compounds"
# PATH = "F:/RITA/TESIS/datos_para_tesis/LA_FE/Metabolomics_Teresa"
# file = "ListMetabolomics_SMILES_canSMILES_notfoundAmano.xlsx"
# col_name = "Compound"
# SEP = ";"
PATH = "G:/Otros ordenadores/Mi_PC_Trabajo/RITA/"
file = "compuestos_Pepe.xlsx"
col_name = "Name"
# SEP = ";"


os.chdir(PATH)

if file.endswith("xlsx"):
    df = pd.read_excel(file)
elif file.endswith("csv"):
    df = pd.read_csv(file, sep=SEP)

names_list = df[col_name].to_list()


print("Obtaining compounds from names")
compounds_list = []
for i, name in enumerate(names_list):
    print(f"Compound {str(i)} of {str(len(names_list))}")
    compounds = pcp.get_compounds(name, "name")
    if len(compounds) >= 1:
        compounds_list.append(compounds[0])
    else:
        compounds_list.append("Not found")

print("Obtaining CIDs and SMILES")
cids_list = []
can_smis_list = []
for i, compound in enumerate(compounds_list):
    print(f"Compound {str(i)} of {str(len(compounds_list))}")
    if compound == "Not found":
        cids_list.append("Not found")
        can_smis_list.append("Not found")
    else:
        cids_list.append(compound.cid)
        can_smis_list.append(compound.canonical_smiles)

df["canon_SMILES"] = can_smis_list    
df["CID"] = cids_list


file_nam_list = file.split(".")
file2 = file_nam_list[0] + "_CID." + file_nam_list[1]

if file.endswith("xlsx"):
    df.to_excel(file2, index=False)
elif file.endswith("csv"):
    df.to_csv(file2, sep=SEP, index=False)

print(f"Table saved in file {file2}")