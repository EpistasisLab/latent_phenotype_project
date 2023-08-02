import numpy as np
import pandas as pd
import os
import torch
import torch.nn as nn
import shap
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from copy import deepcopy as COPY
from tqdm import tqdm
from matplotlib import pyplot as plt
from captum.attr import IntegratedGradients
from sklearn.linear_model import LinearRegression as LR
from sklearn.cluster import KMeans
from umap import UMAP
from time import time
import argparse
import pdb

X_df = pd.read_csv("X.txt", delimiter = "\t", header = 0)
X_cols = np.array(X_df.columns)
X_cols[X_cols == "22001-0.0"] = "is_male"
X_df.columns = X_cols
X = X_df[X_cols[X_cols != "eid"]].to_numpy() 
X_cols = X_cols[X_cols != "eid"]
is_PC = np.array([col[0:2] == "PC" for col in X_cols])
PCs = X[:, is_PC]

prefix = "../step1_get_phenotypes_simple/"
death_ICD_data_all = pd.read_csv(prefix + "death_ICD_data.txt", delimiter = "\t", low_memory = False)
filtered_eids = pd.read_csv("y.txt", delimiter = "\t", usecols = ["eid"])
death_ICD_data_all = death_ICD_data_all.merge(filtered_eids, on = "eid", how = "inner")
data_cols = death_ICD_data_all.columns.to_numpy()[death_ICD_data_all.columns.to_numpy() != "eid"]
death_ICD_data = death_ICD_data_all[data_cols]

ICD_seqs = [row[pd.isnull(row) == False].tolist() for row in death_ICD_data.to_numpy()]
ICD_seqs = [row if len(row) > 0 else [""] for row in ICD_seqs]
ICD_codes = death_ICD_data.to_numpy().reshape(-1)
ICD_codes = ICD_codes[pd.isnull(ICD_codes) == False]
ICD_codes = np.unique(ICD_codes)

ICD_codes2 = np.array(ICD_codes).reshape(-1, 1)
Y = np.zeros((len(ICD_seqs), len(ICD_codes)), dtype = bool)
for i, seq in enumerate(ICD_seqs): Y[i, :] = np.any(np.array(seq).reshape(-1, 1) == np.array(ICD_codes), axis = 0)
out_dim = len(Y[0])

N = 20
main_ICDs = ["I251", "I251", "I259", "I259"]
main_ICDs += ["I259", "I251", "I251", "I251"]
main_ICDs += ["I251", "I251", "I251", "I251"]
main_ICDs += ["I251", "I251", "I251"]

for i in range(15): 

    fname2 = "final_model_shapley_values/step7j_network" + str(i) + "_shapley_value_errors.txt"
    SV_errs = pd.read_csv(fname2, sep = "\t", header = None)
    print("max err for phenotype " + str(i) + ":" + str(np.max(SV_errs.to_numpy())))

    fname_out = "step7k_shapley_plot_pheno" + str(i) + ".png"
    if not os.path.exists(fname_out):
        info_rows = (np.any(Y == 1, axis = 1))
        fname1 = "final_model_shapley_values/step7j_network" + str(i) + "_shapley_values.txt"
        SVs = pd.read_csv(fname1, sep = "\t", header = None)
        MASVs = np.mean(np.abs(SVs.to_numpy())[info_rows], axis = 0)
        sorted_inds = np.flip(np.argsort(MASVs))
        top_codes = ICD_codes[sorted_inds[0:N]]
        top_vals = MASVs[sorted_inds[0:N]]
        top_cols = SVs.to_numpy()[info_rows][:, sorted_inds[0:N]]
        top_features = Y[info_rows][:, sorted_inds[0:N]].astype(float)
        shap.summary_plot(top_cols, top_features, feature_names = top_codes)
        plt.savefig(fname_out)
        plt.clf()

        '''
        pdb.set_trace()

        plt.hist(top_cols[top_features[:, 2] == 1, 2], 100)
        plt.savefig("aaa.png")
        plt.clf()
        main_index = np.where(top_codes == "R13")[0][0]
        main_corrs = [spearmanr(top_cols[:, main_index], vec)[0] for vec in top_cols.T]
        main_corrs[main_index] = 0
        main_corrs = np.sort(main_corrs)[-5:]
        main_corr_ICD = top_codes[np.argmax(main_corrs)]
        '''
        

pdb.set_trace()

'''

#TODO:
#1) analyze SVs 
#2) reconstruct HF phenotype and get residuals
#3) correct for genetic PCs

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

# put fake phenotype generation code here if replacing

phenotypes = pd.DataFrame(phenotypes)
phenotypes["eid"] = X_df["eid"].astype(int)
phenotypes = phenotypes[["eid"] + phenotypes.columns[phenotypes.columns != "eid"].tolist()]
phenotypes.to_csv("logistic_SVD_output/phenotypes" + str(k) + ".txt", sep = "\t", header = True, index = False)
'''

