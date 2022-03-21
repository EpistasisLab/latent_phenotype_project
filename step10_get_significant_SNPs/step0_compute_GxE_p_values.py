import numpy as np
import pandas as pd
import argparse
import os
import fcntl
from functools import reduce
from tqdm import tqdm
from copy import deepcopy as COPY
from scipy.stats import linregress
import statsmodels.api as sm
from scipy.stats import t
from scipy.stats import norm
from bed_reader import open_bed
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from sklearn.linear_model import LinearRegression
import pdb
import matplotlib as mpl
from scipy.stats import norm
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


def perm_test(g0, p0, E0, M):
    
    valued_indices = np.logical_or(np.isnan(p0), np.isnan(g0)) == False
    g = g0[valued_indices]
    g1_vals = (g == 1)
    g2_vals = (g == 2)
    p = p0[valued_indices]
    E1 = E0[valued_indices]
    E = COPY(E1).astype(float)
    ge = np.array([g, E, np.ones(len(g))]).T

    p_vals = []
    if np.all(np.sum(g2_vals, axis = 0) > 1000) and chr != "Y": 
        X = np.array([g1_vals, g2_vals, (g1_vals)*E, (g2_vals)*E, E, np.ones(len(g1_vals))]).T
        for i in tqdm(range(M)):
            np.random.shuffle(E)
            X[:, 2] = (g1_vals)*E
            X[:, 3] = (g2_vals)*E
            X[:, 4] = E
            
            Betas = np.linalg.lstsq(X, p, rcond = None)[0]
            enc1, enc2 = Betas[0] + Betas[2]*E[g1_vals], Betas[1] + Betas[3]*E[g2_vals]  
            g[g1_vals] = enc1     
            g[g2_vals] = enc2  

            ge[:, 0] = g
            ge[:, 1] = E
            p_vals.append(np_OLS(ge, ge[:, 1:], p))

    else:
        g1_vals += g2_vals
        X = np.array([g1_vals, (g1_vals)*E, E, np.ones(len(g1_vals))]).T
        for i in tqdm(range(M)):
            np.random.shuffle(E)
            X[:, 1] = (g1_vals)*E
            X[:, 2] = E

            Betas = np.linalg.lstsq(X, p, rcond = None)[0]
            enc = Betas[0] + Betas[1]*E[g1_vals]  
            g[g1_vals] = enc 

            ge[:, 0] = g
            ge[:, 1] = E
            p_vals.append(np_OLS(ge, ge[:, 1:], p))

    return(p_vals)

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

#  python step0_compute_GxE_p_values.py --chr 4 --pheno_index 7 --rsID rs2466455

env_dict = {'_gender':'22001-0.0', '_smoking':'pack-years', '_alcohol':'annual-consumption', '_exercise':'874-average'}
env_name = env_dict[name]
path = "../step9_regress_phenotypes_against_SNPs/env_factors_cleaned.txt"
env_data = pd.read_csv(path, delimiter = "\t", usecols = ["22001-0.0", env_name])
is_male = env_data["22001-0.0"].to_numpy(dtype = "bool")
env_factor = env_data[env_name].to_numpy()
env_factor[np.isnan(env_factor)] = np.nanmedian(env_factor)

QTL_pheno_path = "../step9_regress_phenotypes_against_SNPs/QTL_phenotypes.txt"
QTL_phenotypes = pd.read_csv(QTL_pheno_path, delimiter = "\t", header = None).to_numpy()

best_hits_path = "hits_QTL" + name + "/QTL_hits_chr" + chr + ".txt" 
best_hits = pd.read_csv(best_hits_path, delimiter = "\t")
row = np.logical_and(best_hits["rsID"] == rsID, best_hits["pheno_index"] == pheno_index)
rsID_info = best_hits[row]

new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
new_fam_paths = np.array([new_fam_path_prefix + i + ".fam" for i in chromosomes])
old_fam_path = "../step4_remove_relatives/UKB_samples_unrelated.fam"
new_fams = [pd.read_csv(path, delim_whitespace = True, header = None)[1].to_numpy() for path in new_fam_paths]
new_fam_intersect = reduce(np.intersect1d, new_fams)

main_path = new_fam_paths[chr_index]
new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
is_in_intersect = np.isin(new_fam_main, new_fam_intersect)
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

file = pd.read_csv("hits_QTL" + name + "/QTL_hits_chr" + chr + ".txt", delimiter = "\t")
cond = np.logical_and(file["rsID"] == rsID, file["pheno_index"] == int(pheno_index))
info = file.loc[cond, "p_main"].values[0]
p_cutoff = -np.log10(info)

p_vals = perm_test(g0, p0, E0, 20000)
W = -np.log10(p_vals)
def T(z): return((p_cutoff  - (A + B*(1/(g))*(np.exp((g)*z[0])-1)*np.exp(h*((z[0])**2)/2)))**2)

sampled_p = []
for i in tqdm(range(10000)):
    linargs = [0.1, 99.9, 1000]
    W2 = np.random.choice(W, len(W))
    A, B, g, h, success = estimate_tukey_params(W2, linargs)
    sampled_p.append(norm.sf(minimize(T, 3).x)[0])

pdb.set_trace()
p_val = np.percentile(sampled_p, 95)
fname = rsID + name + "_" + str(pheno_index) + "_" + chr + "_" + str(p_val) + ".txt"
file = open("hits_GxE_p_vals/" + fname, "w")
file.close()

'''
z = np.random.normal(0, 1, int(100000000))
fitted_TGH = (A + B*(1/(g))*(np.exp((g)*z)-1)*np.exp(h*(z**2)/2))
plt.hist(fitted_TGH, bins = 200, density = True, fc=(0, 0, 1, 0.5))
plt.hist(W, bins = 200, density = True, fc=(1, 0, 0, 0.5))
plt.savefig("thing.png")
plt.clf()
'''