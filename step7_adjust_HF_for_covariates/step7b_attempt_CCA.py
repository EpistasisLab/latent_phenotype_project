import numpy as np
import pandas as pd
import torch
from sklearn.cross_decomposition import CCA
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
for i in range(len(X[0])): X[np.isnan(X[:, i]), i] = np.nanmedian(X[:, i])
y_df = pd.read_csv("y.txt", delimiter = "\t", header = 0)
y_cols = y_df.columns
Y = y_df[y_cols[1:]].to_numpy(dtype = float)

cca = CCA(n_components = 15, max_iter = 2000)
cca.fit(Y, X)
Y_c, X_c = cca.transform(Y, X)
Y_reconstructed = cca.inverse_transform(Y_c)
percent_reconstructed = pearsonr(Y.reshape(-1), Y_reconstructed.reshape(-1))[0]**2
HF_residuals = Y[:, 4] - Y_reconstructed[:, 4]
residual_env_corr_p_vals = [pearsonr(HF_residuals, X[:, i])[1] for i in range(len(X[0]))]
m = np.sum(np.array(residual_env_corr_p_vals) < 0.5)
M = len(residual_env_corr_p_vals)
status1 = binom_test(m, M, 0.5, "greater") > 0.05
status2 = np.min(residual_env_corr_p_vals) > 0.05/len(residual_env_corr_p_vals)
ICD_env_corr_remains = status1 or status2

residuals = []
i = 0
for x, y in zip(X_c.T, Y_c.T):
    i += 1
    m, b, r, p, err = linregress(x,y)
    y_est = m*x + b
    res = y - y_est
    residuals.append(res)
    plt.hist(res, bins = 1000)
    plt.title("feature " + str(i) + " histogram")
    plt.savefig("feature" + str(i) + ".png")
    plt.clf()
    
residuals.append(HF_residuals)
plt.hist(HF_residuals, bins = 1000)
plt.title("HF residuals histogram")
plt.savefig("HF_residuals.png")
Y_c_res = pd.DataFrame(np.array(residuals).T)
Y_c_res["eid"] = X_df["eid"].astype(int)
Y_c_res = Y_c_res[["eid"] + Y_c_res.columns[Y_c_res.columns != "eid"].tolist()]
Y_c_res.to_csv("phenotypes.txt", sep = "\t", header = False, index = False)