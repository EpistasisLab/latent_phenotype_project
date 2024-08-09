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
#------------------------------------------------------------------------------------------------------------------------------------
# start of step 7 methodology change for reviewers
#------------------------------------------------------------------------------------------------------------------------------------
X_df["eid"] = X_df["eid"].astype(int)
path = "/home/greggj/pleiotropy_and_GxE/step1_get_my_phenotypes_complex/reformatted_fields.txt"
townsend_index = pd.read_csv(path, delimiter = "\t", usecols = ["eid", "189-average"])
townsend_index.columns = ["eid", "townsend index"]
X_df = X_df.merge(townsend_index, on = "eid", how = "inner")
X_df.loc[np.isnan(X_df["townsend index"]), "townsend index"] = np.mean(X_df["townsend index"])
X_df.loc[np.isnan(X_df["age"]), "age"] = np.mean(X_df["age"])
#------------------------------------------------------------------------------------------------------------------------------------
# end of step 7 methodology change for reviewers
#------------------------------------------------------------------------------------------------------------------------------------

X_cols = np.array(X_df.columns)
X_cols[X_cols == "22001-0.0"] = "is_male"
X_df.columns = X_cols
X = X_df[X_cols[X_cols != "eid"]].to_numpy() 
X_cols = X_cols[X_cols != "eid"]
is_PC = np.array([col[0:2] == "PC" for col in X_cols])

#------------------------------------------------------------------------------------------------------------------------------------
# start of step 7 methodology change for reviewers
#------------------------------------------------------------------------------------------------------------------------------------
# age, gender, and townsend index features are now included in the regression. 
# also, only PCs 1-6 are included
is_PC2 = np.isin(X_cols, ["age", "is_male", "townsend index"])
is_PC += is_PC2
include_if_PC = np.isin(X_cols, ["PC" + str(i) for i in (np.arange(14) + 7)]) == False
is_PC *= include_if_PC
#------------------------------------------------------------------------------------------------------------------------------------
# end of step 7 methodology change for reviewers
#------------------------------------------------------------------------------------------------------------------------------------
PCs = X[:, is_PC]

y_df = pd.read_csv("y.txt", delimiter = "\t", header = 0)
y_cols = y_df.columns
Y = y_df[y_cols[1:]].to_numpy(dtype = float)

X_fname = "logistic_SVD_output/raw_X" + str(k) + ".txt"
Y_fname = "logistic_SVD_output/raw_Y" + str(k) + ".txt"
if not (os.path.exists(X_fname) and os.path.exists(Y_fname)):

    max_iter = 10000
    tol = 0
    Y_c, Y_rots, Y_rec, scree, runtime, errs, deltas = compute_components(Y, k, max_iter, tol, XY = None)
    metadata = pd.DataFrame(np.array([max_iter*[scree], max_iter*[runtime], errs, deltas], dtype = str).T)
    metadata.columns = ["scree", "runtime", "errs", "deltas"]
    metadata.to_csv("logistic_SVD_output/metadata" + str(k) + ".txt", sep = "\t", header = True, index = False)

    raw_X = pd.DataFrame(Y_c)
    raw_Y = pd.DataFrame(Y_rots)
    raw_X.to_csv(X_fname, sep = "\t", header = False, index = False)
    raw_Y.to_csv(Y_fname, sep = "\t", header = False, index = False)

else: 
    Y_c = pd.read_csv(X_fname, delimiter = "\t", header = None).to_numpy()
    Y_rots = pd.read_csv(Y_fname, delimiter = "\t", header = None).to_numpy()
    XY = np.matmul(Y_c[:, 0:k], Y_rots[0:k, :])
    XY[XY < -88] = -88
    Y_rec = 1/(1 + np.exp(-XY))

any_HF_index = np.where(y_cols[1:].to_numpy()=="any_HF")[0][0]
HF_residuals = Y[:, any_HF_index] - Y_rec[:, any_HF_index]
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
phenotypes.to_csv("logistic_SVD_output/phenotypes" + str(k) + ".txt", sep = "\t", header = True, index = False)