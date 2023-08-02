import numpy as np
import pandas as pd
import os
import statsmodels.api as sm
from sklearn.decomposition import PCA
import pdb

PCs_df = pd.read_csv("X.txt", delimiter = "\t", usecols = ["eid"] + ["PC" + str(i) for i in range(1, 21)])
HF_df = pd.read_csv("y.txt", delimiter = "\t", usecols = ['eid', 'any_HF'])
if np.all(PCs_df["eid"] == HF_df["eid"]) == False:
    print("error: eids from X.txt and Y.txt are not in the same order when they should be")
    exit()

HF_cols = np.array(HF_df.columns)[np.array(HF_df.columns)!= "eid"]
PCs_cols = np.array(PCs_df.columns)[np.array(PCs_df.columns)!= "eid"]
HF, PCs =  HF_df[HF_cols].to_numpy().T, PCs_df[PCs_cols].to_numpy()
PCs2 = np.concatenate((PCs, np.ones((len(PCs), 1))), axis = 1)

betas = []
for y in HF:
    log_reg = sm.Logit((y.astype(int)).reshape(-1,1), PCs2).fit()
    betas.append(log_reg.params)

betas = pd.DataFrame(np.array(betas).T)
betas.columns = HF_cols
betas.to_csv("binary_HF_PC_betas.txt", sep = "\t", header = True, index = False)
HF_df.to_csv("binary_HF_values.txt", sep = "\t", header = True, index = False)
PCs_df.to_csv("binary_HF_genetic_PCs.txt", sep = "\t", header = True, index = False)

