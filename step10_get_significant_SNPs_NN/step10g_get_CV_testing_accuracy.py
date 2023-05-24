import numpy as np
import pandas as pd
import statsmodels.api as sm
import argparse
from itertools import combinations
from functools import reduce
from copy import deepcopy as COPY
from bed_reader import open_bed
from matplotlib import pyplot as plt
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import wilcoxon
from scipy.stats import rankdata
from scipy.stats import mannwhitneyu as mwu
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier as GBC
from sklearn.model_selection import KFold
from scipy.stats import chi2
from scipy.stats import combine_pvalues
from scipy.stats import cauchy
import statsmodels.api as sm
from tqdm import tqdm
import os
import pdb

# TODO
# 2) divide training indices into train and test sets for nested 10-fold CV inner loop
# 3) transform SNPs only over training-training indices
# 4) determine and write into a file optimal hyperparameters for each fold
# 5) compute and document the generalization error as done before

parser = argparse.ArgumentParser()
parser.add_argument('--CV', nargs = 1, type = int, action = "store", dest = "CV")
args = parser.parse_args()
CV = args.CV[0]
CV_inds = pd.read_csv("step10f_train_inds.txt", delimiter = "\t", usecols = [CV])
CV_inds = CV_inds.to_numpy().reshape(-1)
CV_inds = np.sort(CV_inds[np.isnan(CV_inds) == False].astype(int))
holdout_inds = pd.read_csv("step10f_test_inds.txt", delimiter = "\t", usecols = [CV])
holdout_inds = holdout_inds.to_numpy().reshape(-1)
holdout_inds = np.sort(holdout_inds[np.isnan(holdout_inds) == False].astype(int))

# I simply transform the whole dataset solely from the training samples.
def transform_genotypes(genotypes, phenotypes, Env_factor, train_ind, impute = True):
    
    g_new, p, E1 = np.zeros(genotypes.shape), COPY(phenotypes), COPY(Env_factor)
        
    for i in range(len(genotypes[0])):

        g = COPY(genotypes[:, i])
        if impute:
            g[np.isnan(g)] = np.nanmedian(g)
        else:
            val_indices = np.logical_or(np.isnan(g), np.isnan(p)) == False
            g, p, E1 = g[val_indices], p[val_indices], E1[val_indices]

        g2_vals = (g == 2)
        if np.all(np.sum(g2_vals, axis = 0) > 1000) and chr != "Y": 

            g1_vals = (g == 1)
            X = np.array([g1_vals, g2_vals, (g1_vals)*E1, (g2_vals)*E1, E1]).T
            X2 = np.array([g1_vals, g2_vals]).T
            y = (p - np.mean(p))
            
            Betas = LinearRegression().fit(X[train_ind, :], y[train_ind]).coef_
            enc1, enc2 = Betas[0] + Betas[2]*E1[g1_vals], Betas[1] + Betas[3]*E1[g2_vals]  
            g_old = COPY(g)
            g[g1_vals] = enc1     
            g[g2_vals] = enc2   
            ge = np.array([g, E1, np.ones(len(g))]).T

        else:
            g[g2_vals] = 1
            g1_vals = (g == 1)
            X = np.array([g1_vals, (g1_vals)*E1, E1]).T
            y = (p - np.mean(p))

            Betas = LinearRegression().fit(X[train_ind, :], y[train_ind]).coef_
            enc = Betas[0] + Betas[1]*E1[g1_vals]  
            g_old = COPY(g)
            g[g1_vals] = enc 

        if impute:
            g_new[:, i] = g
        else:
            g_new[val_indices, i] = g
            g_new[val_indices == False, i] = np.nan

    return(g_new)

def nanlinregress(x, y):
    xval_indices = (np.isnan(x) == False)
    yval_indices = (np.isnan(y) == False)
    val_indices = np.logical_and(xval_indices, yval_indices)
    x, y = x[val_indices], y[val_indices]
    return(linregress(x, y))

GxE_effects_fnames = [name for name in os.listdir("hits_GxE_p_vals")]
GxE_effects_info = [name[:-4].split("_") for name in GxE_effects_fnames]
chr = [info[-1] for info in GxE_effects_info]
pheno_index = [int(info[-2]) for info in GxE_effects_info]
env_factor = [info[-3] for info in GxE_effects_info]
rsID = ["_".join(info[:-3]) for info in GxE_effects_info]
GxE_effects_pvals = [pd.read_csv("hits_GxE_p_vals/" + name, sep = "\t", header = None) for name in GxE_effects_fnames]
GxE_effects_pvals = [df.to_numpy()[0].tolist() for df in GxE_effects_pvals]
pvals0 = [i[0] for i in GxE_effects_pvals]
pvals1 = [i[1] for i in GxE_effects_pvals]
pvals2 = [i[2] for i in GxE_effects_pvals]
pvals_complete_set = np.array([i for i in zip(pvals0, pvals1)])
GxE_effects_df = pd.DataFrame(zip(rsID, env_factor, pheno_index, chr, pvals0, pvals1, pvals2))
GxE_effects_df.columns = ["rsID", "env_factor", "pheno_index", "chr", "pEDGE2", "pEDGE", "p_joined"]

# formula for corrected m (g*) in https://pubmed.ncbi.nlm.nih.gov/22588983/
# formula for ICC: https://en.wikipedia.org/wiki/Intraclass_correlation
pheno_path = "../step9_regress_phenotypes_against_SNPs_NN/QTL_phenotypes.txt"
X = pd.read_csv(pheno_path, delimiter = "\t", header = None).to_numpy()[:, :15]
N, K = len(X), len(X[0])
xn, x, s2 = np.mean(X, axis = 1), np.mean(X), np.var(X)
ICC = (K*np.sum((xn - x)**2)/(N*(K - 1)*s2)) - (1/(K - 1))
m = (K + 1) - (1 + (K - 1)*ICC)
pb = 5E-8/(m*6)

# getting main effects with classical testing
suffixes = ["smoking", "alcohol", "exercise", "gender"]
prefixes = ["hits_QTL_" + suf for suf in suffixes]
all_paths = []
for prefix in prefixes: all_paths += [prefix + "/" + path for path in os.listdir(prefix)]
env_factors = [path.split("/")[0].split("_")[-1] for path in all_paths]
all_main_effects = [pd.read_csv(path, delimiter = "\t") for path in all_paths]
for i in range(len(all_main_effects)): all_main_effects[i]["env_factor"] = env_factors[i]
all_main_effects = pd.concat(all_main_effects)
main_effects = (all_main_effects[all_main_effects["p_null2"] <= pb])
main_effects = main_effects.drop_duplicates("rsID")

# getting permutation effects
GxE_effects = GxE_effects_df[GxE_effects_df["pEDGE2"] <= pb].sort_values(by = "pEDGE2")
EDGE_effects = GxE_effects_df[GxE_effects_df["pEDGE"] <= pb].sort_values(by = "pEDGE")
mixed_inds = np.logical_and(GxE_effects_df["p_joined"] <= pb, GxE_effects_df["pEDGE"] > pb)
mixed_inds = np.logical_and(mixed_inds, GxE_effects_df["pEDGE2"] > pb)
mixed_effects = GxE_effects_df[mixed_inds].sort_values(by = "p_joined")

eid_path = "../step9_regress_phenotypes_against_SNPs_NN/QTL_phenotypes_eids.txt"
eids = pd.read_csv(eid_path, delimiter = "\t", header = None)
eids.columns = ["eid"]

names = ["_alcohol", "_smoking", "_gender"]
env_IDs = ['annual-consumption', 'pack-years', '22001-0.0']
exercise_IDs = ['874-average', '894-average', '914-average']
path = "../step7_adjust_HF_for_covariates_NN/env_factors_for_step9.txt"
env_data = pd.read_csv(path, delimiter = "\t", usecols = ["eid"] + env_IDs + exercise_IDs)

exercise = env_data[exercise_IDs].to_numpy()
pca = PCA(n_components = 1)
pca.fit(exercise)
exercise_score = pca.transform(exercise).reshape(-1)
env_data["exercise_score"] = exercise_score
env_data = env_data[["eid"] + env_IDs + ["exercise_score"]]
env_data = env_data.merge(eids, on = "eid", how = "inner")
del env_data["eid"]
env_data = env_data.to_numpy()

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
new_fam_paths = np.array([new_fam_path_prefix + i + ".fam" for i in chromosomes])
old_fam_path = "../step4_remove_relatives/UKB_samples_unrelated.fam"
new_fams = [pd.read_csv(path, delim_whitespace = True, header = None)[1].to_numpy() for path in new_fam_paths]
new_fam_intersect = reduce(np.intersect1d, new_fams)

prefixes = ["hits_QTL" + name + "/" for name in names]
path_sets = [[prefix + path for path in os.listdir(prefix)] for prefix in prefixes]
files = [pd.concat([pd.read_csv(path, delimiter = "\t") for path in set]) for set in path_sets]

plink_prefixes = ["significant_SNPs_plink_files" + name + "/" for name in names]
plink_path_sets = [[prefix + path for path in os.listdir(prefix)] for prefix in plink_prefixes]
plink_bed_sets = [[path for path in set if ".bed" in path] for set in plink_path_sets]
plink_bim_sets = [[path for path in set if ".bim" in path] for set in plink_path_sets]
plink_fam_sets = [[path for path in set if ".fam" in path] for set in plink_path_sets]
plink_chr_sets = [np.array([path.split("chr")[1].split(".")[0] for path in set]) for set in plink_bed_sets]

chromosomes_used = pd.concat([file[["chr"]] for file in files])["chr"]
chromosomes_used = np.unique(chromosomes_used.to_numpy(dtype = str))
phenotypes = pd.concat([file[["pheno_index"]] for file in files])["pheno_index"]
phenotypes = np.unique(phenotypes.to_numpy())

if not os.path.isfile("dummy.bim"):
    dummy = pd.DataFrame(2*[6*["dummy"]])
    dummy.to_csv("dummy.bim", sep = "\t", header = False, index = False)
if not os.path.isfile("dummy.bed"):
    file = open("dummy.bed", "w")
    file.close()

#-----------------------------------------------------------------------------------------------------------------------
#
# For inner loop train/test ind generation
# 
#-----------------------------------------------------------------------------------------------------------------------

path = "../step7_adjust_HF_for_covariates_NN/y.txt"
data = pd.read_csv(path, delimiter = "\t") 
data = data.merge(eids, on = "eid", how = "inner")
y = data['any_HF'].to_numpy(dtype = int)
y = y[CV_inds]

N_total = len(CV_inds)
N_cases = np.sum(y)
N_conts = N_total - N_cases

inds_cases = CV_inds[y == 1]
inds_conts = CV_inds[y == 0] 

kf_cases = KFold(n_splits = 10, shuffle = True)
kf_conts = KFold(n_splits = 10, shuffle = True)
train_ind_sets, test_ind_sets = [], []
for N, kf, inds in zip([N_cases, N_conts], [kf_cases, kf_conts], [inds_cases, inds_conts]):
    kf.get_n_splits(np.ones((N, 2)))
    ind_sets = [ind_set for ind_set in kf.split(np.ones((N, 2)))]
    train_ind_sets.append(pd.DataFrame([inds[i[0]] for i in ind_sets]).T)
    test_ind_sets.append(pd.DataFrame([inds[i[1]] for i in ind_sets]).T)

z1 = pd.concat(train_ind_sets).to_numpy()
z2 = pd.concat(test_ind_sets).to_numpy()
check1 = np.array([len(np.unique(z1[np.isnan(z1[:, i]) == False, i])) + len(np.unique(z2[np.isnan(z2[:, i]) == False, i])) for i in range(10)])
check2 = len(np.unique(z1.reshape(-1)[np.isnan(z1.reshape(-1)) == False]))
if not (np.all(check1 == N_total) and check2 == N_total):
    print("there is a bug in the cross validation index division")
    pdb.set_trace()

train_ind_sets = [np.sort(z1[np.isnan(z1[:, i]) == False, i]).astype(int) for i in range(10)]
test_ind_sets = [np.sort(z2[np.isnan(z2[:, i]) == False, i]).astype(int) for i in range(10)]

#-----------------------------------------------------------------------------------------------------------------------
#
# End of inner loop train/test ind generation
# 
#-----------------------------------------------------------------------------------------------------------------------

GBC_r2_val_sets, HP_sets, r2_LR_CV = [], [], []

for train_inds, test_inds in zip(train_ind_sets, test_ind_sets): 
    EDGE2_geno_pheno_sets = [[] for i in range(len(phenotypes))]
    geno_pheno_sets = [[] for i in range(len(phenotypes))]

    # Output for this loop: EDGE2_geno_pheno_sets
    # EDGE2_geno_pheno_sets contains one list per phenotype
    # each such list contains each SNP for main effects and each EDGE2 transformed SNP for GxE effects
    for chr in tqdm(chromosomes_used):

        chr_index = np.where(chr == np.array(chromosomes))[0][0]
        main_path = new_fam_paths[chr_index]
        new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
        is_in_intersect = np.isin(new_fam_main, eids.to_numpy().reshape(-1))
        sorted_main_indices = np.argsort(new_fam_main)
        sorted_indices = sorted_main_indices[is_in_intersect[sorted_main_indices]]

        main_effects_chr = main_effects.loc[main_effects["chr"] == chr, ["rsID", "pheno_index"]]
        prefix = "significant_SNPs_plink_files_gender/"
        main_effects_bim_path = prefix + "significant_SNPs_chr" + chr + ".bim"
    
        # appends every main effect SNP in this loop
        if os.path.isfile(main_effects_bim_path):
            main_effects_rsIDs = pd.read_csv(main_effects_bim_path, delim_whitespace = True, header = None)
            main_effects_bed_path = prefix + "significant_SNPs_chr" + chr + ".bed"
            main_effects_bed_opener = open_bed(main_effects_bed_path, count_A1 = False, num_threads = 1)
            for i in phenotypes:
                main_effects_chr_pheno = main_effects_chr.loc[main_effects_chr["pheno_index"] == i, "rsID"]
                main_effects_ind = np.where(main_effects_rsIDs.isin(main_effects_chr_pheno.to_numpy()))[0]
                main_effects_geno = main_effects_bed_opener.read(np.s_[sorted_indices,  main_effects_ind]).T
                for g in main_effects_geno:
                    g[np.isnan(g)] = np.nanmedian(g)
                    if np.mean(g) > 1: g = 2 - g
                    if np.sum(g == 2) < 1000: g[g == 2] = 1
                    EDGE2_geno_pheno_sets[i].append(g)
                    geno_pheno_sets[i].append(g)
    
        chr_indices = [np.where(set == chr)[0] for set in plink_chr_sets]
        all_E_names, num_E = [name[1:] for name in names], range(len(chr_indices))
        E_name_index_pairs = combinations(num_E, 2)

        QTL_pheno_path = "../step9_regress_phenotypes_against_SNPs_NN/QTL_phenotypes.txt"
        QTL_phenotypes = pd.read_csv(QTL_pheno_path, delimiter = "\t", header = None).to_numpy()
        QTL_phenotypes = (QTL_phenotypes - np.nanmean(QTL_phenotypes, axis = 0))/(np.nanstd(QTL_phenotypes, axis = 0))
    
        chr_bed_paths = [np.array(plink_bed_sets[k])[chr_indices[k]] for k in num_E]
        chr_bed_paths = ['' if len(path) == 0 else path[0] for path in chr_bed_paths]
        chr_bim_paths = [np.array(plink_bim_sets[k])[chr_indices[k]] for k in num_E]
        chr_bim_paths = ["dummy.bim" if len(path) == 0 else path[0] for path in chr_bim_paths]

        chr_bed_openers = ['' if len(path) == 0 else open_bed(path, count_A1 = False, num_threads = 1) for path in chr_bed_paths]
        chr_bim_files = [pd.read_csv(path, delim_whitespace = True, header = None) for path in chr_bim_paths]
        results_partial = GxE_effects[GxE_effects["chr"] == chr]   

        # appends every GxE EDGE2 encoded SNP 
        for i in phenotypes:
            p = COPY(QTL_phenotypes[:, i])
            p2 = COPY(QTL_phenotypes[:, i])
            p[np.isnan(p)] = np.nanmedian(p)
            geno_pheno_sets 
            results =  results_partial[results_partial["pheno_index"] == i]
            E_names = np.unique(results["env_factor"])
            if len(results) > 0:
                hits_subsets = [results.loc[results["env_factor"] == name, "rsID"] for name in all_E_names]
                hits_indices = [np.where(chr_bim_files[k][1].isin(hits_subsets[k].to_numpy()))[0] for k in num_E]
                genotypes = [np.array([]) if len(hits_indices[k]) == 0 else chr_bed_openers[k].read(np.s_[:,  hits_indices[k]]).T for k in num_E]
                for j in num_E:
                    rsIDs_j = hits_subsets[j].to_numpy()
                    geno_sub = genotypes[j]
                    if len(geno_sub) > 0:
                        geno_sub =  geno_sub[:, sorted_indices].T
                        afs = np.nanmean(geno_sub, axis = 0)
                        af_is_not_maf = (afs > 1)
                        geno_sub[:, af_is_not_maf] = 2 - geno_sub[:, af_is_not_maf]

                        E = env_data[:, j]
                        geno_EDGE2 = transform_genotypes(geno_sub, p, E, train_inds)
                        E_perm = np.random.choice(E, len(E), replace = False)
                        geno_EDGE2_null = transform_genotypes(geno_sub, p, E_perm, train_inds)
                        for g1, g2, g2_null, rsID in zip(geno_sub.T, geno_EDGE2.T,  geno_EDGE2_null.T, rsIDs_j):
                            EDGE2_geno_pheno_sets[i].append(g2)
                            geno_pheno_sets[i].append(g1)                        

    # again, easy to implement CV by training model params only on the training indices
    r2_vals = []
    p_vals = []
    PRS_scores = []
    r2_vals_null = []
    p_vals_null = []
    #for i in tqdm(phenotypes[:-1]):
    for i in tqdm(phenotypes):
        N = len(env_data)
        ph = COPY(QTL_phenotypes[:, i])
        ph[np.isnan(ph)] = np.nanmedian(ph)
        if len(EDGE2_geno_pheno_sets[i]) > 0:
            geno_EDGE2 = np.concatenate([np.array(EDGE2_geno_pheno_sets[i]).T, env_data, np.ones((N, 1))], axis = 1)
            geno_null =  np.concatenate([np.array(geno_pheno_sets[i]).T, np.ones((N, 1))], axis = 1)
            for k in range(len(geno_null[0])): geno_null[np.isnan(geno_null[:, k]), k] = np.nanmedian(geno_null[:, k])
            model1 = sm.OLS(ph, geno_EDGE2).fit()
            model2 = sm.OLS(ph[train_inds], geno_EDGE2[train_inds, :]).fit()
            L = len(env_data[0]) + 1
            beta = model2.params[:-L]
            ph_PRS = np.matmul(geno_EDGE2[:, :-L], beta)
            r, p = pearsonr(ph, ph_PRS)
            PRS_scores.append(ph_PRS)
            r2_vals.append(r**2)
            p_vals.append(p)
            model2 = sm.OLS(ph, geno_null).fit()
            r2_vals_null.append(model2.rsquared) 

    path = "../step7_adjust_HF_for_covariates_NN/y.txt"
    data = pd.read_csv(path, delimiter = "\t") 
    data = data.merge(eids, on = "eid", how = "inner")
    y = data['any_HF'].to_numpy(dtype = int)
    X = np.array(PRS_scores).T
    X = (X - np.mean(X, axis = 0))/np.std(X, axis = 0)
    X_train, X_test = X[train_inds], X[test_inds]
    y_train, y_test = y[train_inds], y[test_inds]

    X2, y2 = [], []
    for X, y in zip([X_train, X_test], [y_train, y_test]):
        case_inds = np.where(y == 1)[0]
        control_inds = np.where(y == 0)[0]
        keepers1 = case_inds
        keepers2 = np.random.choice(control_inds, size = len(case_inds), replace = False)
        keepers = np.concatenate([keepers1, keepers2])
        X2.append(X[keepers])
        y2.append(y[keepers])

    LRs = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]
    max_depth = [1, 2, 3, 4, 5, 6, 7]
    param_combos = []
    GBC_r2_vals = []
    for k in tqdm(max_depth):
        for i in LRs:
            param_combos.append((i, k))
            clf = GBC(n_estimators = 100, learning_rate = i, max_depth = k)
            model_info = clf.fit(X2[0], y2[0])
            pps = model_info.predict_proba(X2[1])[:, 1]
            GBC_r2_vals.append(pearsonr(pps, y2[1])[0]**2)

    clf2 = LogisticRegression(C = 1E20, solver = 'sag')
    model_info2 = clf2.fit(X2[0], y2[0])
    pps2 = model_info2.predict_proba(X2[1])[:, 1]
    r2_LR = pearsonr(pps2, y2[1])[0]**2

    GBC_r2_val_sets.append(GBC_r2_vals)
    HP_sets.append(param_combos)
    r2_LR_CV.append(r2_LR)
    
GBC_r2_val_means = np.mean(GBC_r2_val_sets, axis = 0)
best_index = np.argmax(GBC_r2_val_means)
best_params = HP_sets[0][best_index]

y = data['any_HF'].to_numpy(dtype = int)
X = np.array(PRS_scores).T
X = (X - np.mean(X, axis = 0))/np.std(X, axis = 0)
X_train, X_test = X[CV_inds], X[holdout_inds]
y_train, y_test = y[CV_inds], y[holdout_inds]

X3, y3 = [], []
for Xi, yi in zip([X_train, X_test], [y_train, y_test]):
    case_inds = np.where(yi == 1)[0]
    control_inds = np.where(yi == 0)[0]
    keepers1 = case_inds
    keepers2 = np.random.choice(control_inds, size = len(case_inds), replace = False)
    keepers = np.concatenate([keepers1, keepers2])
    X3.append(Xi[keepers])
    y3.append(yi[keepers])

clf = GBC(n_estimators = 100, learning_rate = best_params[0], max_depth = best_params[1])
model_info = clf.fit(X3[0], y3[0])
pps = model_info.predict_proba(X3[1])[:, 1]
GBC_r2 = pearsonr(pps, y3[1])[0]**2

clf2 = LogisticRegression(C = 1E20, solver = 'sag')
model_info2 = clf2.fit(X3[0], y3[0])
pps2 = model_info2.predict_proba(X3[1])[:, 1]
LR_r2 = pearsonr(pps2, y3[1])[0]**2

data = pd.DataFrame(2*[[best_params[0], best_params[1], GBC_r2, LR_r2]])
data.columns = ["learning rate", "max depth", "R squared (GBC nested CV)", "R squared (LR nested CV)"]
fname = "CV_params/fold" + str(CV) + "_nested.txt"
data.to_csv(fname, sep = "\t", header = True, index = False)

LRs = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]
max_depth = [1, 2, 3, 4, 5, 6, 7]
param_combos = []
GBC_r2_vals = []
for k in tqdm(max_depth):
    for i in LRs:
        param_combos.append((i, k))
        clf = GBC(n_estimators = 100, learning_rate = i, max_depth = k)
        model_info = clf.fit(X3[0], y3[0])
        pps = model_info.predict_proba(X3[1])[:, 1]
        GBC_r2_vals.append(pearsonr(pps, y3[1])[0]**2)

data = pd.DataFrame(np.concatenate([param_combos, np.array(GBC_r2_vals).reshape(-1,1)], axis = 1))
data.columns = ["learning rate", "max depth", "R squared (GBC)"]
fname = "CV_params/fold" + str(CV) + ".txt"
data.to_csv(fname, sep = "\t", header = True, index = False)

