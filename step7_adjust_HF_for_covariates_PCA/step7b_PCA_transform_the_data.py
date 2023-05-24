import numpy as np
import pandas as pd
import torch
import os
from bed_reader import open_bed
from sklearn.cross_decomposition import CCA
from sklearn.decomposition import FastICA
from sklearn.decomposition import PCA
from sklearn.decomposition import KernelPCA as kPCA
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
from tqdm import tqdm
import pdb

from step0_autoenc_library import AE
from step0_autoenc_library import AE_penalty
from step0_autoenc_library import k_fold_CV
from step0_autoenc_library import retrain_network

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

def is_field(col_name, fields):

    """
    Purpose
    -------
    Parameters
    ----------
    Returns
    -------
    """

    status = False
    column_has_0th_instance = False

    # These lines assume standard colname notation of "field-X.Y", where X is the instance number, and Y is the rep number.
    partial_col_names = [field + "-" for field in fields]
    colname_has_partial = np.any([col_name[0:(len(part))] == part for part in partial_col_names])
    if colname_has_partial: column_has_0th_instance = col_name.split("-")[1][0] == "0"

    # Also, if there is only 1 instance and 1 rep, then the colname is just "field"
    if column_has_0th_instance or col_name in fields: status = True
    
    return(status)

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

X_cors = np.corrcoef(X.T) 
X_cor_cors = np.corrcoef(X_cors)**2 - np.eye(len(X_cors[0])) 
X_pairs = np.array(np.where(X_cor_cors > 0.5)).T    
X_label_sets_ind = link_pairs(X_pairs)      
X_label_sets = ["_".join(X_cols[ind]) for ind in X_label_sets_ind]
X_cors2 = (np.corrcoef(X_cors)**2)*np.sign(X_cors)

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

pcaX = (PCA(n_components = 15)).fit(X)
pcaY = (PCA(n_components = 15)).fit(Y)
X_c = pcaX.transform(X)
Y_c = pcaY.transform(Y) 
X_rots = (pcaX.components_).T
Y_rots = (pcaY.components_).T

HF_R2s = [0]
all_R2s = [0]
any_HF_R2s = [0]
K_vec = np.arange(16)
for k in tqdm(K_vec[1:]):
    raw_A = Y_c[:, 0:k]
    raw_B = Y_rots[:, 0:k].T
    Y_rec = np.matmul(raw_A, raw_B)
    percent_logPCA_Y = np.array([pearsonr(Y[:, i], Y_rec[:, i])[0]**2 for i in range(len(Y[0]))])
    all_R2s.append(np.mean(percent_logPCA_Y))
    HF_R2s.append(np.mean(percent_logPCA_Y[0:5]))
    any_HF_R2s.append(np.mean(percent_logPCA_Y[4]))

plt.plot(K_vec, HF_R2s, "r-", label = "HF_reconstruction")
plt.plot(K_vec, any_HF_R2s, "g-", label = "any_cause_HF_reconstruction")
plt.plot(K_vec, all_R2s, "c-", label = "ALL_reconstruction")
plt.plot(K_vec, np.ones(len(K_vec)), "k-", label = "uper limit")
plt.legend()
plt.savefig("PCA_reconstruction_plot.png")

X_rots  = X_rots/np.max(np.abs(X_rots), axis = 0) + 1E-9
Y_rots  = Y_rots/np.max(np.abs(Y_rots), axis = 0) + 1E-9

pdb.set_trace()

N = 10
top_phenotypes = []
top_phenotype_scores = []
top_env_factors = []
top_env_factor_scores = []
for i in range(15):
    Y_roti = Y_rots[:, i]
    Y_importance, Y_names = get_group_importance(Y_label_sets_ind, Y_label_sets, Y_cors2, Y_roti, y_cols[1:])
    top_Y_ind = np.flip(np.argsort(np.abs(Y_importance))[-N:])
    top_phenotypes.append(Y_names[top_Y_ind])
    top_phenotype_scores.append(Y_importance[top_Y_ind])

    X_roti = X_rots[:, i]
    X_importance, X_names = get_group_importance(X_label_sets_ind, X_label_sets, X_cors2 , X_roti, X_cols)
    top_X_ind = np.flip(np.argsort(np.abs(X_importance))[-N:])
    top_env_factors.append(X_names[top_X_ind])
    top_env_factor_scores.append(X_importance[top_X_ind])

y_group_names = Y_label_sets + ["ALT_" + l + "_ALT" for l in Y_label_sets]
y_group_names2 = ["heart_fail", "A_Fibril", "heart_ills", "e_varices"]
y_group_names2 += ["heart_fail_het", "A_Fibril_het", "heart_ills_het", "e_varices_het"]
x_group_names = X_label_sets + ["ALT_" + l + "_ALT" for l in X_label_sets]
x_group_names2 = ["insurance", "metabolism", "pollution", "exercise", "abuse", "rape", "20526"]
x_group_names2 += ["20527", "20529", "20530", "BP_or_chol_meds", "vitamins", "eats_added_sugar", "highschool_grad"]
x_group_names2 += [name + "_het" for name in x_group_names2]

phenotypes_att = []
env_factors_att = []
for i in range(15): 
    pheno = top_phenotypes[i]
    for j in range(len(y_group_names)): pheno[pheno == y_group_names[j]] = y_group_names2[j]
    scores = top_phenotype_scores[i].astype('U10')
    #scores[scores == "1.0"] = "1.00000000000000"
    #scores[scores == "-1.0"] = "-1.00000000000000"
    phenotypes_att += [pheno.astype('U20'), top_phenotype_scores[i].astype('U10')]

    env = top_env_factors[i]
    for j in range(len(x_group_names)): env[env == x_group_names[j]] = x_group_names2[j]
    scores = top_env_factor_scores[i].astype('U10')
    #scores[scores == "1.0"] = "1.00000000000000"
    #scores[scores == "-1.0"] = "-1.00000000000000"
    env_factors_att += [env.astype('U20'), top_env_factor_scores[i].astype('U10')]

pd.DataFrame(phenotypes_att).to_csv("top_phenotypes.txt", sep = "\t", header = False, index = False)
pd.DataFrame(env_factors_att).to_csv("top_env_factors.txt", sep = "\t", header = False, index = False)


X_rec = pcaX.inverse_transform(X_c)
Y_rec = pcaY.inverse_transform(Y_c) 

pdb.set_trace()
HF_residuals = Y[:, 4] - Y_rec[:, 4]

percent_PCA_Y = np.array([pearsonr(Y[:, i], Y_rec[:, i])[0]**2 for i in range(len(Y[0]))])
percent_PCA_X = np.array([pearsonr(X[:, i], X_rec[:, i])[0]**2 for i in range(len(X[0]))])

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

# including a fake phenotype solely from 10 SNPs in chr6 (no hits there) to confirm the detection of a known hit. 

prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr6"
col_indices = np.cumsum(40*[18000])
fam_file = pd.read_csv(prefix + ".fam", delim_whitespace = True, header = None)
bim_file = pd.read_csv(prefix + ".bim", delim_whitespace = True, header = None)
genotypes = open_bed(prefix  + ".bed", count_A1 = False, num_threads = 1).read(np.s_[:, col_indices])
non_mafs = np.where(np.nanmean(genotypes, axis = 0)/2 >= 0.5)
for i in non_mafs: genotypes[:, i] = 2 - genotypes[:, i] 
mafs = np.nanmean(genotypes, axis = 0)/2
good_cols = np.logical_and(mafs >= 0.3, np.sum(np.isnan(genotypes), axis = 0) < 5000)
genotypes = genotypes[:, good_cols]
# just shows minimal LD
for i in range(len(genotypes[0])): genotypes[np.isnan(genotypes[:, i]), i] = np.nanmedian(genotypes[:, i]) 
print(np.max(np.corrcoef(genotypes.T)**2 - np.eye(len(genotypes[0]))))
used_col_indices = col_indices[good_cols]
used_rsIDs = bim_file[1].to_numpy()[used_col_indices]

def get_effect(x, y, E):
    X = np.array([x, E, (x*E), np.ones(len(E))]).T
    basic_model_p = sm.OLS(y, X[:, [0, 3]]).fit().pvalues[0]
    model = sm.OLS(y, X)
    model_sub = sm.OLS(y, X[:, [0,1,3]])
    model_results = model.fit()
    model_sub_results = model_sub.fit()
    pval = model_results.compare_lr_test(model_sub_results)[1]
    rsq = model_results.rsquared
    ydev = np.abs(y - np.median(y))
    p_vQTL = pearsonr(ydev, x)[1]
    return([basic_model_p, rsq, pval, p_vQTL])

'''
ratios = []
N = len(genotypes)
P = len(genotypes[0])
is_male = X_df[["eid", "is_male"]]
fam_eids = fam_file[[1]]
fam_eids.columns = ["eid"]
males = fam_eids.merge(is_male, on = "eid", how = "inner")["is_male"].to_numpy(dtype = bool)
females = males == False
for i in tqdm(range(1000)):
    E = np.random.normal(0, 1, (N, P))
    E_reducer = np.zeros(E.shape)
    E_reducer[:, 0:4] = 1
    GxE_reducer = np.zeros(E.shape)
    GxE_reducer[:, 0:2] = 1
    noise = np.random.normal(0,1, (N, P))
    eff = np.sum(0.028*genotypes + 0.022*genotypes*E*GxE_reducer + 0.4*E*E_reducer + 0.45*noise, axis = 1)
    m = get_effect(genotypes[males, 0], eff[males], E[males, 0])[2]
    f = get_effect(genotypes[females, 0], eff[females], E[females, 0])[2]
    ratios.append(np.max([np.log10(m)/np.log10(f), np.log10(f)/np.log10(m)]))
'''

N = len(genotypes)
P = len(genotypes[0])
E = np.random.normal(0, 1, (N, P))
E_reducer = np.zeros(E.shape)
E_reducer[:, 0:4] = 1
GxE_reducer = np.zeros(E.shape)
GxE_reducer[:, 0:2] = 1
noise = np.random.normal(0,1, (N, P))
eff = np.sum(0.028*genotypes + 0.022*genotypes*E*GxE_reducer + 0.4*E*E_reducer + 0.45*noise, axis = 1)
ests = np.array([get_effect(genotypes[:, i], eff, E[:, i]) for i in range(P)])

phenotypes = pd.DataFrame(phenotypes)
phenotypes["eid"] = X_df["eid"].astype(int)
phenotypes = phenotypes[["eid"] + phenotypes.columns[phenotypes.columns != "eid"].tolist()]
fake_pheno = fam_file.loc[:, [1]]
fake_pheno[16] = eff
fake_pheno.columns = ["eid", 16]
phenotypes = phenotypes.merge(fake_pheno, on = "eid", how = "outer")
phenotypes.to_csv("phenotypes.txt", sep = "\t", header = True, index = False)

info = pd.DataFrame(ests.T)
info.columns = used_rsIDs
info.to_csv("simulated_phenotype_info.txt", sep = "\t", header = True, index = False)

fake_env_factor = fam_file.loc[:, [1]]
fake_env_factor[[100, 101, 102, 103, 104, 105]] = E 
fake_env_factor.columns = ["eid", 100, 101, 102, 103, 104, 105]
fake_env_factor.to_csv("fake_env_factors.txt", sep = "\t", header = True, index = False)


