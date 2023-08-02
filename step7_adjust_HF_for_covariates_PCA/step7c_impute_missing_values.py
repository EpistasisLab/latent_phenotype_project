import numpy as np
import pandas as pd
import os
import pdb 
import argparse
from tqdm import tqdm
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import yeojohnson as yj
from copy import deepcopy as COPY

parser = argparse.ArgumentParser()
parser.add_argument('--nn', nargs = 1, type = str, action = "store", dest = "nn")
args = parser.parse_args()
nn = float((args.nn)[0])

Y = pd.read_csv("y.txt", delimiter = "\t")
ICD_coders = np.sum(Y.to_numpy(dtype = int)[:, 1:], axis = 1) > 1
P = pd.read_csv("phenotypes_PCA.txt", delimiter = "\t")
P2 = COPY(P)
for i in np.arange(16).astype(str):
    pheno = (P[i].to_numpy() - np.mean(P[i].to_numpy()))/np.std(P[i].to_numpy())
    pheno2 = yj(pheno)[0]
    P2.loc[:, i] = pheno2
    if nn == 0.1:
        plt.hist(pheno2, bins = 100)
        plt.savefig("phenotype" + i + ".png")
        plt.clf()
        plt.hist(pheno2[ICD_coders], bins = 100)
        plt.savefig("phenotype" + i + "_with_ICD_codes.png")
        plt.clf()

X_all = pd.read_csv("X.txt", delimiter = "\t")
# These are features that weren't considered for anything before
E_general = pd.read_csv("../step1_get_phenotypes_simple/phenotype_info.txt", delimiter = "\t", low_memory = False)
E_general = E_general.merge(Y[["eid"]], on = "eid", how = "inner")
E_cols = E_general.columns.to_numpy()

# Generating complete IC code data
E_cols = E_general.columns.to_numpy()
colfields = np.array([name.split("-")[0] for name in E_cols])
ICD_codes_1 = E_general.loc[:, E_cols[np.isin(colfields, ["eid", "40002"])]]
ICD_codes_2 = E_general.loc[:, E_cols[np.isin(colfields, ["eid", "40001"])]]
ICD_codes_3 = E_general.loc[:, E_cols[np.isin(colfields, ["eid", "41270"])]]
ICD_codes_all = ICD_codes_1.merge(ICD_codes_2, how = "inner", on = "eid")
ICD_codes_all = ICD_codes_all.merge(ICD_codes_3, how = "inner", on = "eid")
del ICD_codes_all["eid"]
ICD_codes_all = ICD_codes_all.to_numpy(dtype = str).reshape(-1)
ICD_codes_all = np.unique(ICD_codes_all[ICD_codes_all != "nan"])

if not os.path.exists("ICD_codes_bin.txt"):
    ICD_codes_bin = E_general.loc[:, E_cols[np.isin(colfields, ["40001", "40002", "41270"])]].to_numpy(dtype = str)
    ICD_codes_bin = np.array([np.isin(ICD_codes_all, row) for row in ICD_codes_bin])
    ICD_codes_bin = pd.DataFrame(ICD_codes_bin)
    ICD_codes_bin.columns = ICD_codes_all
    ICD_codes_bin["eid"] = E_general["eid"]
    ICD_codes_bin = ICD_codes_bin[["eid"] + ICD_codes_all.tolist()]
    # too many values without this step 
    val_counts = np.sum(ICD_codes_bin.to_numpy(), axis = 0)
    ICD_codes_bin = ICD_codes_bin[ICD_codes_bin.columns[val_counts > 30]]
    (ICD_codes_bin.astype(int)).to_csv("ICD_codes_bin.txt", sep = "\t", header = True, index = False)
else:
    ICD_codes_bin = pd.read_csv("ICD_codes_bin.txt", delimiter = "\t")

# removing ICD codes, self reported illnesses, and related info
E_cols = [name for name in E_cols if "41280" not in name]
E_cols = [name for name in E_cols if "41270" not in name]
E_cols = [name for name in E_cols if "40002" not in name]
E_cols = [name for name in E_cols if "40001" not in name]
E_cols = [name for name in E_cols if "20002" not in name]
# removing complicated consult reason info
E_cols = [name for name in E_cols if "41245" not in name]
E_cols = [name for name in E_cols if "41246" not in name]
E_cols = [name for name in E_cols if "41249" not in name]
# removing cols that we already have (some are under different names)
modified_cols = ['34-0.0', '52-0.0', '129-0.0', '130-0.0', "699-0.0", "738-0.0"]
modified_cols += ["21000-0.0", "21003-0.0", "22006-0.0"]
X_cols = X_all.columns.to_numpy()[X_all.columns.to_numpy() != "eid"]
E_cols = [name for name in E_cols if name not in modified_cols]
E_cols = [name for name in E_cols if name not in X_cols]
E_general = E_general[E_cols]
E_general[E_general.astype(float) < 0] = np.nan
X_all = X_all.merge(E_general, on = "eid", how = "inner")

cols_to_impute = ['eid', 'pack-years', 'annual-consumption', '874-average', '894-average', '914-average']
# some exercise values are imputed differently.
# all exercise values with all three missing are imputed 0 because it likely means they reported not exercising (31007/54627 missings). 
# all moderate and heavy exercise values with both missing are imputed 0 (79530/97162 misings)
# if moderate or heavy exercise is nonzero, light walking is likely also nonzero
# if heavy exercise is nonzero, moderate exercise is also likely nonzero
other_cols = [name for name in X_all.columns if "PC" not in name]
other_cols = np.setdiff1d(other_cols, cols_to_impute[1:])
other_cols = other_cols[other_cols != "C2"]
X = X_all.loc[:, cols_to_impute]
X_other = COPY(X_all.loc[:, other_cols])
X_other.loc[X_other["C1"] == -1, "C1"] = np.nan

unique_val_sets = [np.unique(X_all.loc[:, col]) for col in other_cols]
unique_val_counts = np.array([len(S[np.isnan(S) == False]) for S in unique_val_sets])
#cols_to_log_transform = other_cols[unique_val_counts > 10]
#cols_to_log_transform = cols_to_log_transform[cols_to_log_transform != "eid"]
#for col in cols_to_log_transform:
#    X_other.loc[:, col] = np.log(X_other[col] - np.nanmin(X_other[col]) + np.nanmean(X_other[col]))
#    plt.hist(X_other.loc[np.isnan(X_other[col]) == False, col].to_numpy(), bins = 100)
#    plt.savefig("info_features/info_feature_" + col + ".png")
#    plt.clf()
is_imp_col_names = [col + "_is_imp" for col in other_cols]
is_imp_cols = pd.DataFrame(np.array([np.isnan(X_other[col].to_numpy()).astype(int) for col in other_cols]).T)
is_imp_cols.columns = is_imp_col_names 
is_imp_cols["eid"] = X_other["eid"]
for col in other_cols: X_other.loc[np.isnan(X_other[col]), col] = np.nanmean(X_other[col].to_numpy())
X_other = X_other.merge(is_imp_cols, on = "eid", how = "inner")
X_other = X_other.merge(ICD_codes_bin, on = "eid", how = "inner")
other_cols = np.setdiff1d(X_other.columns, cols_to_impute[1:])

corr_sets = []
for name in ['pack-years', 'annual-consumption', '874-average', '894-average', '914-average']:
    val_inds = np.isnan(X[name].to_numpy()) == False
    corr_sets.append([pearsonr(X.loc[val_inds, name], X_other.loc[val_inds, col])[0] for col in other_cols])
best_corrs = np.max(corr_sets, axis = 0)
# sometimes "eid" may have a low correlation, and if selected, we don't want it to appear twice
best_cols = np.union1d(["eid"], other_cols[np.abs(best_corrs) > nn])
X_other = X_other.loc[:, best_cols]

no_light_exercise = X_other['864-0.0'] == 0
no_moderate_exercise = X_other['884-0.0'] == 0
no_heavy_exercise = X_other['904-0.0'] == 0
X.loc[no_light_exercise, '874-average'] = 0
X.loc[no_moderate_exercise , '894-average'] = 0
X.loc[no_heavy_exercise , '914-average'] = 0
X_test = COPY(X.dropna())
# log_trans_cols = ["pack-years", "annual-consumption", "874-average", "894-average", "914-average"]
# for col in log_trans_cols: X_test.loc[:, col] = np.log(X_test[col].to_numpy() - np.min(X_test[col].to_numpy()) + np.mean(X_test[col].to_numpy()))

X_test_std = X_test.merge(X_other, on = "eid", how = "inner")
X_test_std_cols = X_test_std.columns.to_numpy()
X_test_std_cols = X_test_std_cols[np.std(X_test_std.to_numpy(), axis = 0) > 0]
X_test_std_cols = X_test_std_cols[X_test_std_cols != "eid"]
X_test_std = X_test_std[X_test_std_cols].to_numpy()
X_test_std = (X_test_std - np.mean(X_test_std, axis = 0))/np.std(X_test_std, axis = 0)
pca = PCA(n_components = 5)
PCs = pca.fit_transform(X_test_std)
nan_indices1 = np.array([np.random.choice(len(X_test), 20000, replace = False) for i in range(len(cols_to_impute) - 1)]).T
nan_indices2 = np.array([np.argsort(X_test[col])[-20000:] for col in cols_to_impute[1:]]).T
nan_indices3 = np.array([np.argsort(X_test[col])[:20000] for col in cols_to_impute[1:]]).T
nan_indices4 = np.array([np.argsort(PCs[:, 0])[-20000:] for col in cols_to_impute[1:]]).T
nan_indices5 = np.array([np.argsort(PCs[:, 0])[:20000] for col in cols_to_impute[1:]]).T
nan_indices6 = np.array([np.argsort(PCs[:, 1])[-20000:] for col in cols_to_impute[1:]]).T
nan_indices7 = np.array([np.argsort(PCs[:, 1])[:20000] for col in cols_to_impute[1:]]).T
nan_indices8 = np.array([np.argsort(PCs[:, 2])[-20000:] for col in cols_to_impute[1:]]).T
nan_indices9 = np.array([np.argsort(PCs[:, 2])[:20000] for col in cols_to_impute[1:]]).T
nan_indices10 = np.array([np.argsort(PCs[:, 3])[-20000:] for col in cols_to_impute[1:]]).T
nan_indices11 = np.array([np.argsort(PCs[:, 3])[:20000] for col in cols_to_impute[1:]]).T
nan_indices12 = np.array([np.argsort(PCs[:, 4])[-20000:] for col in cols_to_impute[1:]]).T
nan_indices13 = np.array([np.argsort(PCs[:, 4])[:20000] for col in cols_to_impute[1:]]).T


info_sets = []

ind_sets = [nan_indices1, nan_indices2, nan_indices3, nan_indices4, nan_indices5]
test_names = ["random", "highest", "lowest", "highest_PC1", "lowest_PC1"]
for inds, name in zip(ind_sets, test_names):
    D = X_test.merge(P2, on = "eid", how = "inner")
    D = D.merge(X_other, on = "eid", how = "inner")
    D_names = D.columns.to_numpy()[D.columns.to_numpy() != "eid"]
    D2 = D[D_names].to_numpy()
    unique_counts = np.array([len(np.unique(col)) for col in D2.T])
    mu = np.nanmean(D2[:, unique_counts > 2], axis = 0)
    sig = np.nanstd(D2[:, unique_counts > 2], axis = 0)
    D2[:, unique_counts > 2] = (D2[:, unique_counts > 2] - mu)/sig
    D3 = COPY(D2)
    for i, set in enumerate(inds.T): D3[set, i] = np.nan
    imputer = IterativeImputer(max_iter=100, random_state=0)
    imputer.fit(D3)
    D_imp = imputer.transform(D3)

    info = [name]
    for i in range(5):
        indices = inds[:, i]
        imputed_vals = np.max([D_imp[indices, i], np.min(D2[:, i])*np.ones(len(indices))], axis = 0)
        real_vals = D2[indices, i]
        imputed_diffs_sum = np.sum(np.abs(real_vals - imputed_vals))
        null_diffs_sum = np.sum(np.abs((real_vals - np.mean(D2[:, i]))))
        var_explained = 1 - (imputed_diffs_sum/null_diffs_sum)
        print(var_explained)
        info.append(var_explained)
    print("DONE")
    info_sets.append(info)

info_sets = pd.DataFrame(info_sets)
info_sets.columns = ["test_name"] + D_names[:5].tolist()
fname = "corr_sets_" + str(nn) + "neighbors.txt"
info_sets.to_csv(fname, sep = "\t", header = True, index = False)