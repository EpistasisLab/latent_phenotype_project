import numpy as np
import pandas as pd
import argparse
import os
import fcntl
from functools import reduce
from copy import deepcopy as COPY
from scipy.stats import linregress
import statsmodels.api as sm
from itertools import combinations as combos
from scipy.stats import t
from scipy.stats import norm
from bed_reader import open_bed
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from scipy.stats import yeojohnson as yj
import pdb
import matplotlib as mpl
from scipy.stats import norm
from tqdm import tqdm
from scipy.stats import chi2
mpl.rcParams['agg.path.chunksize'] = 10000
from statsmodels.stats.outliers_influence import summary_table

def estimate_tukey_params(W, linargs):

    """
    Purpose
    -------
    To estimate tukey parameters for the distribution of deviation scores
    
    
    Parameters
    ----------
    W: the distribution of deviation scores 

    Returns
    -------
    A: tukey location parameter
    B: tukey scale parameter
    g: tukey skew parameter
    h: tukey tail heaviness parameter
    """
    Q_vec = np.percentile(W, [10, 25, 50, 75, 90])
    if len(np.unique(Q_vec)) != 5:
        pdb.set_trace()
        Q_vec = approximate_quantiles(W, [10, 25, 50, 75, 90])

    # rough A estimation
    A = Q_vec[2]

    # rough B estimation
    IQR = Q_vec[3] - Q_vec[1]
    QR2 = Q_vec[4] - Q_vec[0]
    if IQR == 0 or QR2 == 0:
        if len(W[W <  Q_vec[0]]) != 0 and len(W[W > Q_vec[4]]) != 0:
            IQR = np.min(W[W > Q_vec[3]]) - np.max(W[W <  Q_vec[1]])
            QR2 = np.min(W[W > Q_vec[4]]) - np.max(W[W <  Q_vec[0]])
        else:
            IQR = 1.35*np.std(W)
            QR2 = 2.56*np.std(W)
    SK = (Q_vec[4] + Q_vec[0] - 2*Q_vec[2])/QR2
    T = QR2/IQR
    phi = 0.6817766 + 0.0534282*SK + 0.1794771*T - 0.0059595*(T**2)
    B = (0.7413*IQR)/phi

    # rough g estimation
    zv = norm.ppf(0.9, 0, 1)
    UHS = Q_vec[4] - Q_vec[2]
    LHS = Q_vec[2] - Q_vec[0]
    if UHS == 0 or LHS == 0:
        if len(W[W <  Q_vec[0]]) != 0 and len(W[W > Q_vec[4]]) != 0:
            UHS = np.min(W[W >= Q_vec[4]]) - np.max(W[W <  Q_vec[2]])
            LHS = np.min(W[W > Q_vec[2]]) - np.max(W[W <=  Q_vec[0]])
        else:
            UHS = np.max(W) - np.mean(W)
            LHS = np.mean(W) - np.min(W)
            zv = norm.ppf(1 - ((0.1)**np.log10(len(W))), 0, 1)
    g = (1/zv)*np.log(UHS/LHS)

    # rough h estimation
    y = (W - A)/B        
    Q_vec2 = np.percentile(y, [10, 25, 50, 75, 90])
    if len(np.unique(Q_vec2)) != 5:
        pdb.set_trace()
        Q_vec2 = approximate_quantiles(y, [10, 25, 50, 75, 90])
        
    Q_ratio = (Q_vec2[4]*Q_vec2[0])/(Q_vec2[4] + Q_vec2[0])
    if Q_ratio <= 0:
        h = 0
    else:
        h = (2/(zv**2))*np.log(-g*Q_ratio)
    if np.isnan(h) or np.abs(h) == np.inf:
        h = 0

    qi = np.linspace(*linargs)
    pi_real = np.percentile(W, qi)
    zi = norm(0,1).ppf(qi/100)
    old_err = (pi_real - (A + (B/g)*(np.exp(g*zi) - 1)*np.exp(h*(zi**2)/2)))**2
    def tukey_loss(theta):
        A, B, g, h = theta
        pi_est = A + (B/g)*(np.exp(g*zi) - 1)*np.exp(h*(zi**2)/2)
        return(np.mean((pi_est - pi_real)**2))
    theta0 = [A, B, g, h]
    theta_est_data = minimize(tukey_loss, theta0)
    A, B, g, h = theta_est_data.x
    g += 1E-6
    new_err = (pi_real - (A + (B/g)*(np.exp(g*zi) - 1)*np.exp(h*(zi**2)/2)))**2 
    return(A, B, g, h, theta_est_data.success)

def nanlinregress(x, y):
    xval_indices = (np.isnan(x) == False)
    yval_indices = (np.isnan(y) == False)
    val_indices = np.logical_and(xval_indices, yval_indices)
    x, y = x[val_indices], y[val_indices]
    return(linregress(x, y))

# taken from https://online.stat.psu.edu/stat462/node/131/
# and https://online.stat.psu.edu/stat462/node/137/
def np_OLS(X_full, X_less, y):
  
    Betas2_full, SSE_full = np.linalg.lstsq(X_full, y, rcond = None)[0:2]
    Betas2_less, SSE_less = np.linalg.lstsq(X_less, y, rcond = None)[0:2]
    d_full = len(y) - 3
    T = ((SSE_less - SSE_full)/(SSE_full/d_full))**(1/2)
    p_val1 = 2*t(d_full).sf(T)[0]

    '''
    # just confirms my p values are computed correctly
    model = sm.OLS(y, X_full)
    p_val2 = model.fit().pvalues[0]
    if not np.round(p_val1, 6) == np.round(p_val2, 6):
        pdb.set_trace()
    else:
        print(p_val1)
        print(p_val2)
        print("\n\n")
    '''
    return(p_val1)

def compute_EDGE2_p_value(g, p, E, is_male_input):

    g0, p0, E0, is_male = COPY(g), COPY(p), COPY(E), COPY(is_male_input)
    if chr == "Y":
        g0, p0, E0, is_male = g0[is_male], p0[is_male], E0[is_male], is_male[is_male]

    val_inds = np.isnan(g0) == False
    g0, p0, E0, is_male = g0[val_inds], p0[val_inds], E0[val_inds], is_male[val_inds]
    env_is_not_gender = np.any(E0 != is_male)

    if chr == "X" and env_is_not_gender:
        g_set = [g0[is_male], g0[is_male == False]]
        p_set = [p0[is_male], p0[is_male == False]]
        E_set = [E0[is_male], E0[is_male == False]]
        all_vals = zip(g_set, p_set, E_set, [0, 1])
        num_groups = 2
    else:
        all_vals = zip([g0], [p0], [E0], [0])
        num_groups = 1

    EDGE2_p_val_set = []
    for geno, pheno, env, ind in all_vals:

        g2_vals = (geno == 2)
        is_standard = np.sum(g2_vals, axis = 0) > 1000 and chr != "Y"
        if is_standard: 

            g1_vals = (geno == 1)
            X = np.array([g1_vals, g2_vals, (g1_vals)*env, (g2_vals)*env, env, np.ones(len(g1_vals))]).T
            y = (pheno - np.mean(pheno))

            Betas = np.linalg.lstsq(X, y, rcond = None)[0][:-1]
            enc1, enc2 = Betas[0] + Betas[2]*env[g1_vals], Betas[1] + Betas[3]*env[g2_vals]
            geno[g1_vals] = enc1     
            geno[g2_vals] = enc2   

            ge = [geno, env]
            ge = np.array(ge + [np.ones(len(geno))]).T
            EDGE2_p_val = np_OLS(ge, ge[:, 1:], pheno)
            if num_groups == 2:
                EDGE2_p_val_set.append(EDGE2_p_val)

        elif not (chr == "Y" and len(np.unique(E0)) == 1):
            geno[g2_vals] = 1
            g1_vals = (geno == 1)
            X = np.array([g1_vals, (g1_vals)*env, env, np.ones(len(g1_vals))]).T
            y = (pheno - np.mean(pheno))

            Betas = np.linalg.lstsq(X, y, rcond = None)[0][:-1]
            enc = Betas[0] + Betas[1]*env[g1_vals]   
            geno[g1_vals] = enc 

            ge = np.array([geno, env, np.ones(len(geno))]).T
            EDGE2_p_val = np_OLS(ge, ge[:, 1:], pheno)
            if num_groups == 2:
                EDGE2_p_val_set.append(EDGE2_p_val)

    if num_groups == 2:
        EDGE2_p_val = chi2(4).sf(-2*np.sum(np.log(EDGE2_p_val_set), axis = 0))
    return(EDGE2_p_val)

def perm_test(g0, p0, E0, is_male, M):

    E = COPY(E0).astype(float)
    P = COPY(p0).astype(float)
    EDGE2_p_vals = []
    EDGE_p_vals = []
    combined_p_vals = []
    for i in tqdm(range(M)):
            np.random.shuffle(E)
            np.random.shuffle(P)
            EDGE2_p_vals.append(compute_EDGE2_p_value(g0, p0, E, is_male))
    return(EDGE2_p_vals)

# python step0_compute_SV_p_values.py --chr 1 --pheno 12 --rsID rs6698949 --model PCA --count 4
# python step0_compute_SV_p_values.py --chr 1 --pheno 13 --rsID rs11799887 --model PCA --count 5
# python step0_compute_SV_p_values.py --chr 14 --pheno 12 --rsID rs17486830 --model PCA --count 3
parser = argparse.ArgumentParser()
parser.add_argument('--chr', nargs = 1, type = str, action = "store", dest = "chr")
parser.add_argument('--pheno_index', nargs = 1, type = int, action = "store", dest = "pheno_index")
parser.add_argument('--rsID', nargs = 1, type = str, action = "store", dest = "rsID")
parser.add_argument('--model', nargs = 1, type = str, action = "store", dest = "model")
parser.add_argument('--count', nargs = 1, type = int, action = "store", dest = "count")

args = parser.parse_args()
chr = args.chr[0]
pheno_index = args.pheno_index[0]
rsID = args.rsID[0]
model = (args.model)[0]
count = args.count[0]

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
chromosomes = np.array(chromosomes)
chr_index = np.where(chr == chromosomes)[0][0]

if not os.path.exists("hits_GxE_p_vals"):
    os.mkdir("hits_GxE_p_vals")

QTL_pheno_path = "../step7_adjust_HF_for_covariates_PCA/final_model_shapley_values"
QTL_pheno_path += "/NN_phenotype" + str(pheno_index) + "_shapley_values.txt.gz"
QTL_eids_path = "../step9_regress_phenotypes_against_SNPs_NN/QTL_phenotypes_eids.txt"
QTL_eids = pd.read_csv(QTL_eids_path, delimiter = "\t", header = None)
QTL_eids.columns = ["eid"]

SVs = pd.read_csv(QTL_pheno_path, delimiter = "\t", compression = 'gzip')
SVs = SVs.merge(QTL_eids, on = "eid", how = "inner")

X_df = pd.read_csv("X.txt", delimiter = "\t", header = 0)
X_df = X_df.merge(QTL_eids, on = "eid", how = "inner")
X_cols = np.array(X_df.columns)
X_cols = X_cols[X_cols != "eid"]
X = X_df[X_cols].to_numpy() 
is_PC = np.array([col[0:2] == "PC" for col in X_cols])
PCs = X[:, is_PC]

PCs2 = np.concatenate((PCs, np.ones((len(PCs), 1))), axis = 1)
PCsT_PCs = np.matmul(PCs2.T, PCs2)
correction = np.eye(len(PCsT_PCs))*np.min(PCsT_PCs)*1E-6
PCsT_PCs_inv = np.linalg.inv(PCsT_PCs + correction)
coef = np.matmul(PCsT_PCs_inv, PCs2.T)

codes = SVs.columns[SVs.columns.to_numpy() != "eid"]
importance_sums = np.sum(np.abs(SVs.to_numpy()[:, :-1]), axis = 0)
cutoff = np.argmax(np.cumsum(np.flip(np.sort(importance_sums)))/np.sum(importance_sums) > 0.85)
top_codes = codes[np.argsort(importance_sums)][-cutoff:]

'''
from scipy.stats import pearsonr
Z = pd.read_csv("phenotypes_for_step9.txt", delimiter = "\t")
Z = Z.merge(QTL_eids, on = "eid", how = "inner")
Z = Z["12"].to_numpy()
pearsonr(Z, QTL_phenotypes)

vals, counts = np.unique(QTL_phenotypes, return_counts = True)
main_val = vals[np.argmax(counts)]
plt.hist(QTL_phenotypes[QTL_phenotypes != main_val], bins = 200)
plt.savefig("aaa.png")
plt.clf()
'''

env_name = 'pack-years'
path = "../step7_adjust_HF_for_covariates_NN/env_factors_for_step9.txt"
gender_path = "../step7_adjust_HF_for_covariates_NN/X.txt"
gender_data = pd.read_csv(gender_path, delimiter = "\t", usecols = ["eid", "22001-0.0"], dtype = int)
env_data = pd.read_csv(path, delimiter = "\t", usecols = ["eid", env_name])

env_data = env_data.merge(QTL_eids, on = "eid", how = "inner")
gender_data = gender_data.merge(QTL_eids, on = "eid", how = "inner")
is_male = gender_data["22001-0.0"].to_numpy(dtype = "bool")
env_factor = env_data[env_name].to_numpy()

new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
main_path = new_fam_path_prefix + chromosomes[chr_index] + ".fam"
new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
is_in_intersect = np.isin(new_fam_main, QTL_eids.to_numpy().reshape(-1))
sorted_main_indices = np.argsort(new_fam_main)
sorted_indices = sorted_main_indices[is_in_intersect[sorted_main_indices]]

path = "../step10_get_significant_SNPs_" + model + "/significant_SNPs_plink_files_smoking/significant_SNPs_chr" + chr
rsIDs = pd.read_csv(path  + ".bim", delim_whitespace = True, header = None)[1]
rsID_ind =  np.where(rsIDs == rsID)[0][0]
g0_getter = open_bed(path  + ".bed", count_A1 = False, num_threads = 1)
g0 = g0_getter.read(np.s_[sorted_indices,  rsID_ind]).reshape(-1)
if np.nanmean(g0) > 1: g0 = 2 - g0
E0 = env_factor

code_subsets = [i for i in combos(top_codes, count)]
phenotypes, p_EDGE2_vals, p_normal_vals = [], [], []
for subset in tqdm(code_subsets):
    # Y = np.sum(SVs.to_numpy()[:, :-1], axis = 1)
    Y = np.sum(SVs[list(subset)].to_numpy(), axis = 1)
    weights = np.matmul(coef, Y.reshape(-1,1))
    Y_est = np.matmul(PCs2, weights)
    Y_res = Y - Y_est.reshape(-1)
    QTL_phenotypes = yj((Y_res - np.mean(Y_res))/np.std(Y_res))[0]
    phenotypes.append(QTL_phenotypes)
    p_EDGE2_vals.append(-np.log10(compute_EDGE2_p_value(g0, QTL_phenotypes, E0, is_male)))
    p_normal_vals.append(-np.log10(nanlinregress(g0, QTL_phenotypes)[3]))

pdiffs = np.array(p_EDGE2_vals) - np.array(p_normal_vals)
best_ind = np.argmax(pdiffs)
p_cutoff  = p_EDGE2_vals[best_ind]
p_normal = p_normal_vals[best_ind]
p0 = phenotypes[best_ind]
best_ICD_codes = code_subsets[best_ind]

reps = 20000
p_vals = perm_test(g0, p0, E0, is_male, reps)

W = -np.log10(p_vals)
def T(z): return((p_cutoff  - (A + B*(1/(g))*(np.exp((g)*z[0])-1)*np.exp(h*((z[0])**2)/2)))**2)

sampled_p = []
for i in tqdm(range(10000)):
    linargs = [0.1, 99.9, 1000]
    W2 = np.random.choice(W, len(W))
    A, B, g, h, success = estimate_tukey_params(W2, linargs)
    sampled_p.append(norm.sf(minimize(T, 3).x)[0])

p_val = np.percentile(sampled_p, 95)
fname = rsID + "_smoking_" + str(pheno_index) + "_" + chr + "_count_" + str(count) + "_minuslog10p_" + str(-np.log10(p_val))[:5] + ".txt"
top20_ICD10_code_sets = np.array(code_subsets)[np.argsort(pdiffs)[-20:]]
file = pd.DataFrame(top20_ICD10_code_sets)
file["p value proxy"] = np.sort(pdiffs)[-20:]
file = file[["p value proxy"] + list(range(len(top20_ICD10_code_sets[0])))]
file = file.sort_values("p value proxy", ascending = False)
file.to_csv("hits_GxE_p_vals/" + fname, sep = "\t", header = True, index = False)

'''
z = np.random.normal(0, 1, int(100000000))
fitted_TGH = (A + B*(1/(g))*(np.exp((g)*z)-1)*np.exp(h*(z**2)/2))
plt.hist(fitted_TGH, bins = 200, density = True, fc=(0, 0, 1, 0.5))
plt.hist(W, bins = 200, density = True, fc=(1, 0, 0, 0.5))
plt.savefig("thing.png")
plt.clf()
'''