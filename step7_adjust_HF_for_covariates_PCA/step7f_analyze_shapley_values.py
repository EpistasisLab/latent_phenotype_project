import numpy as np
import pandas as pd
import os
import torch
import torch.nn as nn
import shap
import seaborn as sns
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from copy import deepcopy as COPY
from tqdm import tqdm
from matplotlib import pyplot as plt
from captum.attr import IntegratedGradients
from sklearn.linear_model import LinearRegression as LR
from scipy.stats import yeojohnson as yj
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage
from umap import UMAP
from time import time
import argparse
import pdb

Y_df = pd.read_csv("y.txt", delimiter = "\t", header = 0)
all_cols = Y_df.columns.to_numpy()
ICD_cols = all_cols[all_cols != "eid"]
Y = Y_df[ICD_cols].to_numpy(dtype = int)
out_dim = len(Y[0])

if not os.path.isdir("SV_figures"):
    os.mkdir("SV_figures")
N = 30
for i in range(15): 

    fname2 = "final_model_shapley_values/NN_phenotype" + str(i) + "_shapley_value_errors.txt.gz"
    SV_errs = pd.read_csv(fname2, sep = "\t", compression = 'gzip')["error"]
    print("max err for phenotype " + str(i) + ": " + str(np.max(SV_errs.to_numpy())))

    fname_out = "step7f_shapley_plot_pheno" + str(i) + ".png"
    fname_out2 = "step7f_shapley_corr_pheno" + str(i) + ".png"

    info_rows = (np.any(Y == 1, axis = 1))
    fname1 = "final_model_shapley_values/NN_phenotype" + str(i) + "_shapley_values.txt.gz"
    SVs = pd.read_csv(fname1, sep = "\t", compression = 'gzip', low_memory = False)
    eids = SVs["eid"]
    del SVs["eid"]
    SVs.columns = ICD_cols

    MASVs = np.mean(np.abs(SVs.to_numpy())[info_rows], axis = 0)
    sorted_inds = np.flip(np.argsort(MASVs))
    important_codes = ICD_cols[sorted_inds[0:N]]

    all_corrs = np.corrcoef(SVs.to_numpy().T)
    Nr = len(all_corrs)
    cluster_rows = np.any(np.abs(all_corrs - np.eye(len(all_corrs))) >= 0.5, axis = 0)
    cluster_codes = ICD_cols[cluster_rows]

    top_codes = np.union1d(important_codes, cluster_codes)
    SV_subset = SVs.loc[info_rows, top_codes]
    new_cols = np.array([code.split("counts_")[-1] for code in SV_subset.columns])
    new_cols[new_cols == "any_HF"] = "AHF"
    new_cols = [code if len(code) < 4 else code[0:3] + "." + code[3:] for code in new_cols]
    SV_subset.columns = new_cols
    corr = SV_subset.corr(method='pearson')
    abs_dist = lambda u, v: np.sqrt(((np.abs(u)-np.abs(v))**2).sum())
    # abs_dist = lambda u, v: np.sqrt(np.sum((u - v)**2))
    sns.set(font_scale = 3*30/len(corr))
    L = linkage(corr, metric=abs_dist, method = "average")
    # L = linkage(corr, metric = "euclidean", method = "ward")
    M = np.ones((len(corr), len(corr)), dtype = bool)
    # "orientation":"horizontal"
    D = {"aspect":1, "ticks": [-1, 0, 1], "label":"feature shapley values\n correlation coefficient"}
    F = (len(corr), len(corr))
    cg = sns.clustermap(corr, row_linkage = L, col_linkage = L, annot = False, vmin = -1, vmax = 1, figsize = F, cmap = 'RdBu', cbar_kws = D)
    cg.cbar_pos = (0.02, 0.02, 0.02, 0.02)
    cg.ax_row_dendrogram.set_visible(False)
    cg.savefig("SV_figures/" + fname_out2) 
    cg.fig.clf()

    '''
    top_vals = MASVs[sorted_inds[0:N]]
    top_cols = SVs.to_numpy()[info_rows][:, sorted_inds[0:N]]
    SV_subset = SVs.loc[info_rows, top_codes]
    top_features = Y[info_rows][:, sorted_inds[0:N]].astype(float)
    shap.summary_plot(top_cols, top_features, feature_names = top_codes)
    plt.savefig("SV_figures/" + fname_out)
    plt.clf()
    '''
