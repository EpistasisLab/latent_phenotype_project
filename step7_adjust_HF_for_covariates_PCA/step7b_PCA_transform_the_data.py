import numpy as np
import pandas as pd
import torch
import os
from sklearn.decomposition import PCA
from copy import deepcopy as COPY
import pdb

X_df = pd.read_csv("X.txt", delimiter = "\t", header = 0)
X_cols = np.array(X_df.columns)
X_cols[X_cols == "22001-0.0"] = "is_male"
X_df.columns = X_cols
X = X_df[X_cols[X_cols != "eid"]].to_numpy() 
X_cols = X_cols[X_cols != "eid"]
is_PC = np.array([col[0:2] == "PC" for col in X_cols])
PCs = X[:, is_PC]
PCs_df = X_df[["eid"] + ["PC" + str(i) for i in range(1, 21)]]
X =  X[:, is_PC == False]
X_cols = X_cols[is_PC == False]
for i in range(len(X[0])): X[np.isnan(X[:, i]), i] = np.nanmedian(X[:, i])
y_df = pd.read_csv("y.txt", delimiter = "\t", header = 0)
y_cols = y_df.columns
Y = y_df[y_cols[1:]].to_numpy(dtype = float)

pcaY = (PCA(n_components = 15)).fit(Y)
Y_c = pcaY.transform(Y) 
Y_rots = (pcaY.components_).T
Y_rec = pcaY.inverse_transform(Y_c) 

HF_residuals = Y[:, 4] - Y_rec[:, 4]
Y_c_res = np.concatenate([Y_c, HF_residuals.reshape(-1, 1)], axis = 1)
PCs2 = np.concatenate((PCs, np.ones((len(PCs), 1))), axis = 1)
PCsT_PCs = np.matmul(PCs2.T, PCs2)
correction = np.eye(len(PCsT_PCs))*np.min(PCsT_PCs)*1E-6
PCsT_PCs_inv = np.linalg.inv(PCsT_PCs + correction)
coef = np.matmul(PCsT_PCs_inv, PCs2.T)
weights = np.matmul(coef, Y_c_res)
Y_est = np.matmul(PCs2, weights)
phenotypes = Y_c_res - Y_est

phenotypes = pd.DataFrame(phenotypes)
phenotypes["eid"] = X_df["eid"].astype(int)
phenotypes = phenotypes[["eid"] + phenotypes.columns[phenotypes.columns != "eid"].tolist()]
phenotypes.to_csv("phenotypes_PCA.txt", sep = "\t", header = True, index = False)