import numpy as np
import pandas as pd
import os
import argparse
from copy import deepcopy as COPY
import statsmodels.api as sm
from scipy.linalg import sqrtm
from scipy.stats import pearsonr
from scipy.stats import yeojohnson as yj
from scipy.stats import linregress
from matplotlib import pyplot as plt
from scipy.stats import binom_test
from scipy.stats import chi2
from scipy.stats import rankdata
from sklearn.decomposition import TruncatedSVD as SVD
from tqdm import tqdm
from time import time
import pdb

if not os.path.exists("logistic_SVD_output"):
    try: 
        os.mkdir("logistic_SVD_output")
    except:
        pass

parser = argparse.ArgumentParser()
parser.add_argument('--k', nargs = 1, type = str, action = "store", dest = "k")
args = parser.parse_args()
k = int((args.k)[0])

def compute_components(P, k, max_iter, tol, XY = None):
    start = time()
    P = P.astype(np.float32)
    if XY == None:
        svd_output = SVD(n_components = k).fit(P)
        X = svd_output.transform(P)
        Y = svd_output.components_
        XY = np.matmul(X, Y)
        XY[XY < -88] = -88
    num_iter, old_err, err, delta = 0, 1, 0, 1
    errs, deltas = [], []
    for i in tqdm(range(max_iter)):
        G = (1/(1 + np.exp(-XY))) - P
        H = (XY - 4*G)
        svd_output = SVD(n_components = k).fit(H)
        X = svd_output.transform(H)
        Y = svd_output.components_
        XY = np.matmul(X, Y)
        XY[XY < -88] = -88
        err = np.sum((H - XY)**2)/np.sum(H**2)
        delta = old_err - err
        errs.append(err)
        deltas.append(delta)
        old_err = COPY(err)
        num_iter += 1
        if delta > tol:
            continue
    end = time()
    runtime = end - start
    XY = np.matmul(X[:, 0:k], Y[0:k, :])
    XY[XY < -88] = -88
    P_rec = 1/(1 + np.exp(-XY))
    P_null = np.zeros(P.shape) + np.mean(P, axis = 0)
    dev_alt = np.sum(P*np.log(P_rec + 1E-30)) + np.sum((1 - P)*np.log(1 - P_rec + 1E-30))
    dev_null = np.sum(P*np.log(P_null + 1E-30)) + np.sum((1 - P)*np.log(1 - P_null + 1E-30))
    scree_val = 1 - (dev_alt/dev_null)
    return(X, Y, P_rec, scree_val, runtime, errs, deltas)

def link_pairs(pairs):
    vals = np.unique(pairs)
    groups = []
    group_lengths = np.array([])
    for val in vals: 
        subset = np.unique(pairs[pairs[:, 0] == val, :])
        unions = [np.union1d(subset, g) for g in groups]
        union_lengths = np.array([len(u) for u in unions])
        good_union = union_lengths < (group_lengths + len(subset))
        if np.all(union_lengths == (group_lengths + len(subset))):
            groups.append(subset)
            group_lengths = np.array([len(g) for g in groups])
        else:
            for i in range(len(groups)):
                if good_union[i]:
                    groups[i] = unions[i]
            group_lengths = np.array([len(g) for g in groups])
    return(groups)

X_df = pd.read_csv("X.txt", delimiter = "\t", header = 0)
X_cols = np.array(X_df.columns)
X_cols[X_cols == "22001-0.0"] = "is_male"
X_df.columns = X_cols
X = X_df[X_cols[X_cols != "eid"]].to_numpy() 
X_cols = X_cols[X_cols != "eid"]
is_PC = np.array([col[0:2] == "PC" for col in X_cols])
PCs = X[:, is_PC]

y_df = pd.read_csv("y.txt", delimiter = "\t", header = 0)
y_cols = y_df.columns
Y = y_df[y_cols[1:]].to_numpy(dtype = float)

Y_cors = np.corrcoef(Y.T) 
Y_cor_cors = np.corrcoef(Y_cors)**2 - np.eye(len(Y_cors[0])) 
Y_pairs = np.array(np.where(Y_cor_cors > 0.5)).T    
Y_label_sets_ind = link_pairs(Y_pairs)      
Y_label_sets = ["-".join([name.split("s_")[-1] for name in y_cols[1:][ind]]) for ind in Y_label_sets_ind]
Y_cors2 = (np.corrcoef(Y_cors)**2)*np.sign(Y_cors)

def get_group_importance(g_ind, g_labels, cor_cors, rots, labels):

    in_group_ind = np.unique(np.concatenate(g_ind))
    out_group_ind = np.setdiff1d(np.arange(len(rots)), in_group_ind)
    out_group, out_labels = rots[out_group_ind], labels[out_group_ind].tolist()

    group_importance1 = []
    group_importance2 = []
    for ind in g_ind:
        weights1 = np.mean(cor_cors[ind][:, ind], axis = 0)
        weights2 = np.mean(1 - np.abs(cor_cors[ind][:, ind]), axis = 0)
        group_importance1.append(np.sum(weights1*rots[ind]))
        group_importance2.append(np.sum(weights2*np.abs(rots[ind])))

    importance = np.concatenate([out_group, group_importance1, group_importance2])
    alt_g_labels = ["ALT_" + l + "_ALT" for l in g_labels]
    last_labels = np.array(out_labels + g_labels + alt_g_labels)
    return(importance, last_labels)

# pcaY = (PCA(n_components = 15)).fit(Y)
# Y_c = pcaY.transform(Y) 
# Y_rots = (pcaY.components_).T
# Y_rots  = Y_rots/np.max(np.abs(Y_rots), axis = 0) + 1E-9

max_iter = 10000
tol = 0
Y_c, Y_rots, Y_rec, scree, runtime, errs, deltas = compute_components(Y, k, max_iter, tol, XY = None)
metadata = pd.DataFrame(np.array([max_iter*[scree], max_iter*[runtime], errs, deltas], dtype = str).T)
metadata.columns = ["scree", "runtime", "errs", "deltas"]
metadata.to_csv("logistic_SVD_output/metadata" + str(k) + ".txt", sep = "\t", header = True, index = False)

raw_X = pd.DataFrame(Y_c)
raw_Y = pd.DataFrame(Y_rots)
raw_X.to_csv("logistic_SVD_output/raw_X" + str(k) + ".txt", sep = "\t", header = False, index = False)
raw_Y.to_csv("logistic_SVD_output/raw_Y" + str(k) + ".txt", sep = "\t", header = False, index = False)

N = 10
top_phenotypes = []
top_phenotype_scores = []
for i in range(k):
    Y_roti = Y_rots[i, :]
    Y_importance, Y_names = get_group_importance(Y_label_sets_ind, Y_label_sets, Y_cors2, Y_roti, y_cols[1:])
    top_Y_ind = np.flip(np.argsort(np.abs(Y_importance))[-N:])
    top_phenotypes.append(Y_names[top_Y_ind])
    top_phenotype_scores.append(Y_importance[top_Y_ind])

y_group_names = Y_label_sets + ["ALT_" + l + "_ALT" for l in Y_label_sets]
y_group_names2 = ["heart_fail", "A_Fibril", "heart_ills", "e_varices"]
y_group_names2 += ["heart_fail_het", "A_Fibril_het", "heart_ills_het", "e_varices_het"]

phenotypes_att = []
env_factors_att = []
for i in range(k): 
    pheno = top_phenotypes[i]
    for j in range(len(y_group_names)): pheno[pheno == y_group_names[j]] = y_group_names2[j]
    scores = top_phenotype_scores[i].astype('U10')
    phenotypes_att += [pheno.astype('U20'), top_phenotype_scores[i].astype('U10')]

pd.DataFrame(phenotypes_att).to_csv("logistic_SVD_output/top_phenotypes" + str(k) + ".txt", sep = "\t", header = False, index = False)

pdb.set_trace()
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


