import numpy as np
import pandas as pd
import torch
from sklearn.cross_decomposition import CCA
from sklearn.decomposition import PCA
from copy import deepcopy as COPY
from scipy.stats import pearsonr
from scipy.stats import yeojohnson as yj
from scipy.stats import linregress
from matplotlib import pyplot as plt
from scipy.stats import binom_test
import pdb

X_df = pd.read_csv("X.txt", delimiter = "\t", header = 0)
X_cols = np.array(X_df.columns)
X_cols[X_cols == "22001-0.0"] = "is_male"
X_df.columns = X_cols
X = X_df[X_cols[X_cols != "eid"]].to_numpy() 
X_cols = X_cols[X_cols != "eid"]
is_PC = np.array([col[0:2] == "PC" for col in X_cols])
PCs = X[:, is_PC]
X =  X[:, is_PC == False]
X_cols = X_cols[is_PC == False]
for i in range(len(X[0])): X[np.isnan(X[:, i]), i] = np.nanmedian(X[:, i])
y_df = pd.read_csv("y.txt", delimiter = "\t", header = 0)
y_cols = y_df.columns
Y = y_df[y_cols[1:]].to_numpy(dtype = float)

cca = CCA(n_components = 15, max_iter = 2000)
cca.fit(Y, X)
Y_c, X_c = cca.transform(Y, X)

# I swapped their positions by accident, not that it matters except for the notation. 
Y_rots = cca.x_rotations_
X_rots = cca.y_rotations_
N = 10
top_phenotypes = []
top_phenotype_scores = []
top_env_factors = []
top_env_factor_scores = []
for i in range(15):
    top_Y_ind = np.flip(np.argsort(np.abs(Y_rots[:, i]))[-N:])
    top_X_ind = np.flip(np.argsort(np.abs(X_rots[:, i]))[-N:])
    top_phenotypes.append(y_cols[1:][top_Y_ind])
    top_phenotype_scores.append(Y_rots[:, i][top_Y_ind])
    top_env_factors.append(X_cols[top_X_ind])
    top_env_factor_scores.append(X_rots[:, i][top_X_ind])

phenotypes_att = []
for i in range(15): phenotypes_att += [top_phenotypes[i].astype('U15'), top_phenotype_scores[i].astype('U15')]
env_factors_att = []
for i in range(15): env_factors_att += [top_env_factors[i].astype('U15'), top_env_factor_scores[i].astype('U15')]
pd.DataFrame(phenotypes_att).to_csv("top_phenotypes.txt", sep = "\t", header = False, index = False)
pd.DataFrame(env_factors_att).to_csv("top_env_factors.txt", sep = "\t", header = False, index = False)

# this is equivalent to scikit-learn's reconstruction code. 
x_mean = X.mean(axis=0)
x_std = X.std(axis=0, ddof=1)
x_std[x_std == 0.0] = 1.0
y_mean = Y.mean(axis=0)
y_std = Y.std(axis=0, ddof=1)
y_std[y_std == 0.0] = 1.0
Y_rec  = np.matmul(Y_c, cca.x_loadings_.T)*y_std + y_mean
X_rec  = np.matmul(X_c, cca.y_loadings_.T)*x_std + x_mean

# checking to ensure that the HF residuals are uncorrelated to all environmental factors
HF_residuals = Y[:, 4] - Y_rec[:, 4]
residual_env_corr_p_vals = [pearsonr(HF_residuals, X[:, i])[1] for i in range(len(X[0]))]
m = np.sum(np.array(residual_env_corr_p_vals) < 0.5)
M = len(residual_env_corr_p_vals)
status1 = binom_test(m, M, 0.5, "greater") < 0.05
status2 = np.min(residual_env_corr_p_vals) < 0.05/len(residual_env_corr_p_vals)
ICD_env_corr_remains = status1 or status2

percent_CCA_Y = np.array([pearsonr(Y[:, i], Y_rec[:, i])[0]**2 for i in range(len(Y[0]))])
percent_CCA_X = np.array([pearsonr(X[:, i], X_rec[:, i])[0]**2 for i in range(len(X[0]))])

X_residuals = X - X_rec
X_s = (X_residuals - np.mean(X_residuals, axis = 0))/np.std(X_residuals, axis = 0)
pca = PCA(n_components = 15)
pca.fit(X_s)
X_res_pca = pca.transform(X_s)
X_est = pca.inverse_transform(X_res_pca)
percent_PCA_X = np.array([pearsonr(X_s[:, i], X_est[:, i])[0]**2 for i in range(len(X_s[0]))])
percent_rec = percent_CCA_X + (1-percent_CCA_X)*percent_PCA_X

'''
# Manually including some features that look important but aren't accounted for in the CCA or PCA
poorly_explained_features = X_cols[percent_rec < 0.25]
# features included manually:
# age
# income
# 914-average (hours of heavy exercise per week)
'''

additional_features = X_residuals[:, np.isin(X_cols, ["age", "income", "914-average"])]
env_factors = pd.DataFrame(np.concatenate([X_c, X_res_pca, additional_features], axis = 1))
env_factors["eid"] = X_df["eid"].astype(int)
env_factors = env_factors[["eid"] + env_factors.columns[env_factors.columns != "eid"].tolist()]
env_factors.to_csv("env_factors.txt", sep = "\t", header = False, index = False)

Y_c_res = np.concatenate([Y_c, HF_residuals.reshape(-1, 1)], axis = 1)
PCs2 = np.concatenate((PCs, np.ones((len(PCs), 1))), axis = 1)
PCsT_PCs = np.matmul(PCs2.T, PCs2)
correction = np.eye(len(PCsT_PCs))*np.min(PCsT_PCs)*1E-6
PCsT_PCs_inv = np.linalg.inv(PCsT_PCs + correction)
coef = np.matmul(PCsT_PCs_inv, PCs2.T)
weights = np.matmul(coef, Y_c_res)
Y_est = np.matmul(PCs2, weights)
phenotypes = Y_c_res - Y_est

i = 0
for p in phenotypes.T:
    i += 1
    plt.hist(p, bins = 1000)
    plt.title("feature " + str(i) + " histogram")
    plt.savefig("feature" + str(i) + ".png")
    plt.clf()

phenotypes = pd.DataFrame(phenotypes)
phenotypes["eid"] = X_df["eid"].astype(int)
phenotypes = phenotypes[["eid"] + phenotypes.columns[phenotypes.columns != "eid"].tolist()]
phenotypes.to_csv("phenotypes.txt", sep = "\t", header = True, index = False)