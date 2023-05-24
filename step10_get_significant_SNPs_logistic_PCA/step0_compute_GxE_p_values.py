import numpy as np
import pandas as pd
import argparse
import os
import fcntl
from functools import reduce
from copy import deepcopy as COPY
from scipy.stats import linregress
import statsmodels.api as sm
from scipy.stats import t
from scipy.stats import norm
from bed_reader import open_bed
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
import pdb
import matplotlib as mpl
from scipy.stats import norm
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

def compute_EDGE_p_value(g, p, is_male_input, nonlinear = False):

    g0, p0, is_male = COPY(g), COPY(p), COPY(is_male_input)
    env_is_not_gender = np.any(E0 != is_male)
    if chr == "Y":
        g0, p0, is_male = g0[is_male], p0[is_male], is_male[is_male]

    val_inds = np.isnan(g0) == False
    g0, p0, is_male = g0[val_inds], p0[val_inds], is_male[val_inds]

    if chr == "X" and env_is_not_gender:
        g_set = [g0[is_male], g0[is_male == False]]
        p_set = [p0[is_male], p0[is_male == False]]
        all_vals = zip(g_set, p_set, [0, 1])
        num_groups = 2
    else:
        all_vals = zip([g0], [p0], [0])
        num_groups = 1

    EDGE_p_val_set = []
    for geno, pheno, ind in all_vals:

        g2_vals = (geno == 2)
        g1_vals = (geno == 1)
        X = [g1_vals, g2_vals]
        if chr in ["X", "XY"] and env_is_not_gender == False:
            X += [is_male]
        X = np.array(X + [np.ones(len(g1_vals))]).T
        y = (pheno - np.mean(pheno))

        Betas = np.linalg.lstsq(X, y, rcond = None)[0][:-1]
        geno_old = COPY(geno)
        geno[g1_vals] = Betas[0]  
        geno[g2_vals] = Betas[1]

        EDGE_p_val = linregress(geno, pheno)[3]
        if num_groups == 2:
            EDGE_p_val_set.append(EDGE_p_val)

    if num_groups == 2:
        EDGE_p_val = chi2(4).sf(-2*np.sum(np.log(EDGE_p_val_set), axis = 0))
    return(EDGE_p_val)


def perm_test(g0, p0, E0, is_male, status, M):

    E = COPY(E0).astype(float)
    P = COPY(p0).astype(float)
    EDGE2_p_vals = []
    EDGE_p_vals = []
    combined_p_vals = []
    for i in range(M):
            np.random.shuffle(E)
            np.random.shuffle(P)
            EDGE2_p_vals.append(compute_EDGE2_p_value(g0, p0, E, is_male))
            combined_p_vals.append(compute_EDGE2_p_value(g0, P, E, is_male))
            if status != -1: 
                np.random.shuffle(P)
                EDGE_p_vals.append(compute_EDGE_p_value(g0, P, is_male))
    return([EDGE2_p_vals,  EDGE_p_vals, combined_p_vals])


parser = argparse.ArgumentParser()
parser.add_argument('--chr', nargs = 1, type = str, action = "store", dest = "chr")
parser.add_argument('--pheno_index', nargs = 1, type = int, action = "store", dest = "pheno_index")
parser.add_argument('--rsID', nargs = 1, type = str, action = "store", dest = "rsID")
parser.add_argument('--name', nargs = 1, type = str, action = "store", dest = "name")

args = parser.parse_args()
chr = args.chr[0]
pheno_index = args.pheno_index[0]
rsID = args.rsID[0]
name = (args.name)[0]

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
chromosomes = np.array(chromosomes)
chr_index = np.where(chr == chromosomes)[0][0]

if not os.path.exists("hits_GxE_p_vals"):
    os.mkdir("hits_GxE_p_vals")

QTL_pheno_path = "../step9_regress_phenotypes_against_SNPs_logistic_PCA/QTL_phenotypes.txt"
QTL_eids_path = "../step9_regress_phenotypes_against_SNPs_logistic_PCA/QTL_phenotypes_eids.txt"
QTL_phenotypes = pd.read_csv(QTL_pheno_path, delimiter = "\t", header = None).to_numpy()
QTL_eids = pd.read_csv(QTL_eids_path, delimiter = "\t", header = None)
QTL_eids.columns = ["eid"]

#  python step0_compute_GxE_p_values.py --chr 4 --pheno_index 7 --rsID rs2466455

#                                                                   ex             ex               ex          BMI         metabolism   pol 2010
# maybe start with these: ['pack-years', 'annual-consumption', '874-average',  '894-average', '914-average', '21001-0.0', '23105-0.0', '24003-0.0', 'age']

env_dict = {'_gender':'22001-0.0', '_smoking':'pack-years', '_alcohol':'annual-consumption', '_exercise':['874-average', '894-average', '914-average']}
env_name = env_dict[name]
path = "../step7_adjust_HF_for_covariates_logistic_PCA/env_factors_for_step9.txt"
gender_path = "../step7_adjust_HF_for_covariates_logistic_PCA/X.txt"
gender_data = pd.read_csv(gender_path, delimiter = "\t", usecols = ["eid", "22001-0.0"], dtype = int)
if env_name == ['874-average', '894-average', '914-average']:
    env_data = pd.read_csv(path, delimiter = "\t", usecols = (["eid"] + env_name))
    exercise = env_data[env_name].to_numpy()
    pca = PCA(n_components = 1)
    pca.fit(exercise)
    exercise_score = pca.transform(exercise).reshape(-1)
    env_data["exercise_score"] = exercise_score
    env_data = env_data[["eid", "exercise_score"]]
    env_name = "exercise_score"
else:
    env_data = pd.read_csv(path, delimiter = "\t", usecols = ["eid", env_name])

env_data = env_data.merge(QTL_eids, on = "eid", how = "inner")
gender_data = gender_data.merge(QTL_eids, on = "eid", how = "inner")
is_male = gender_data["22001-0.0"].to_numpy(dtype = "bool")
env_factor = env_data[env_name].to_numpy()

best_hits_path = "hits_QTL" + name + "/QTL_hits_chr" + chr + ".txt" 
best_hits = pd.read_csv(best_hits_path, delimiter = "\t")
row = np.logical_and(best_hits["rsID"] == rsID, best_hits["pheno_index"] == pheno_index)
rsID_info = best_hits[row]
p_EDGE2, p_EDGE, p_linear = rsID_info[['p_main', 'p_null1', 'p_null2']].to_numpy().reshape(-1)

new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
main_path = new_fam_path_prefix + chromosomes[chr_index] + ".fam"
new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
is_in_intersect = np.isin(new_fam_main, QTL_eids.to_numpy().reshape(-1))
sorted_main_indices = np.argsort(new_fam_main)
sorted_indices = sorted_main_indices[is_in_intersect[sorted_main_indices]]

path = "significant_SNPs_plink_files" + name + "/significant_SNPs_chr" + chr
rsIDs = pd.read_csv(path  + ".bim", delim_whitespace = True, header = None)[1]
rsID_ind =  np.where(rsIDs == rsID)[0][0]
g0_getter = open_bed(path  + ".bed", count_A1 = False, num_threads = 1)
g0 = g0_getter.read(np.s_[sorted_indices,  rsID_ind]).reshape(-1)
if np.nanmean(g0) > 1: g0 = 2 - g0
p0 = QTL_phenotypes[:, pheno_index]
E0 = env_factor

# -----------------------------------------------------------------------------------------------------
#
# TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Find out why EDGE p value was 0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# python step0_compute_GxE_p_values.py --rsID rs532462830 --pheno_index 8 --chr X --name _smoking
#
# -----------------------------------------------------------------------------------------------------

# pdb.set_trace()
# compute_EDGE_p_value(g0, p0, is_male)
# compute_EDGE2_p_value(g0, p0, E0, is_male)

file = pd.read_csv("hits_QTL" + name + "/QTL_hits_chr" + chr + ".txt", delimiter = "\t")
cond = np.logical_and(file["rsID"] == rsID, file["pheno_index"] == int(pheno_index))
info = file.loc[cond, "p_main"].values[0]
p_cutoff_EDGE2 = -np.log10(p_EDGE2)
if p_EDGE != -1:
    p_cutoff_EDGE = -np.log10(p_EDGE)
else:
    p_cutoff_EDGE = 1

status = p_EDGE
reps = 20000
p_EDGE2_dist, p_EDGE_dist, p_combined_dist = perm_test(g0, p0, E0, is_male, status, reps)

pvals_complete_set = []
for p_vals, p_cutoff in zip([p_EDGE2_dist, p_EDGE_dist, p_combined_dist], [p_cutoff_EDGE2, p_cutoff_EDGE, p_cutoff_EDGE2]):
    if len(p_vals) == reps:
        W = -np.log10(p_vals)
        def T(z): return((p_cutoff  - (A + B*(1/(g))*(np.exp((g)*z[0])-1)*np.exp(h*((z[0])**2)/2)))**2)

        sampled_p = []
        for i in range(10000):
            linargs = [0.1, 99.9, 1000]
            W2 = np.random.choice(W, len(W))
            A, B, g, h, success = estimate_tukey_params(W2, linargs)
            sampled_p.append(norm.sf(minimize(T, 3).x)[0])

        p_val = np.percentile(sampled_p, 95)
        pvals_complete_set.append(p_val)

if len(pvals_complete_set) == 2: pvals_complete_set = [pvals_complete_set[0]] + [p_linear] + [pvals_complete_set[1]]
fname = rsID + name + "_" + str(pheno_index) + "_" + chr + ".txt"
file = pd.DataFrame([pvals_complete_set])
file.to_csv("hits_GxE_p_vals/" + fname, sep = "\t", header = False, index = False)

'''
z = np.random.normal(0, 1, int(100000000))
fitted_TGH = (A + B*(1/(g))*(np.exp((g)*z)-1)*np.exp(h*(z**2)/2))
plt.hist(fitted_TGH, bins = 200, density = True, fc=(0, 0, 1, 0.5))
plt.hist(W, bins = 200, density = True, fc=(1, 0, 0, 0.5))
plt.savefig("thing.png")
plt.clf()
'''