import numpy as np
import pandas as pd
from copy import deepcopy as COPY
from scipy.stats import linregress
from scipy.stats import t
from scipy.stats import norm
from scipy.optimize import minimize
import pdb
import matplotlib as mpl
import argparse
from scipy.stats import chi2
from tqdm import tqdm
mpl.rcParams['agg.path.chunksize'] = 10000

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
    
    chr == "1"
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

    G = COPY(g0).astype(float)
    p_EDGE2 = compute_EDGE2_p_value(g0, p0, E, is_male)
    p_cutoff_EDGE2 = -np.log10(p_EDGE2)
    EDGE2_p_vals = []
    for i in range(M):
            np.random.shuffle(G)
            EDGE2_p_vals.append(compute_EDGE2_p_value(G, p0, E0, is_male))
    return([EDGE2_p_vals, p_cutoff_EDGE2])

# M is the number of permutation test reps (20000)
# B is the number of distribution fitting bootstrap reps (10000)
def perm_test_bootstrap(g0, p0, E0, is_male, M, B):

    p_EDGE2_dist, p_cutoff_EDGE2 = perm_test(g0, p0, E0, is_male, M)
    W = -np.log10(p_EDGE2_dist)
    def T(z): return((p_cutoff_EDGE2  - (A + B*(1/(g))*(np.exp((g)*z[0])-1)*np.exp(h*((z[0])**2)/2)))**2)

    sampled_p = []
    for i in range(B):
        linargs = [0.1, 99.9, 1000]
        W2 = np.random.choice(W, len(W))
        A, B, g, h, success = estimate_tukey_params(W2, linargs)
        sampled_p.append(norm.sf(minimize(T, 3).x)[0])

    rigorous_p_val = np.percentile(sampled_p, 95)
    return(rigorous_p_val)

parser = argparse.ArgumentParser()
parser.add_argument('--i', nargs = 1, type = int, action = "store", dest = "i")
parser.add_argument('--j', nargs = 1, type = int, action = "store", dest = "j")
args = parser.parse_args()
i = args.i[0]
j = args.j[0]

SNPs = pd.read_csv("step12a_SNPs.txt", delimiter = "\t", header = None, usecols = [i])
SNP = SNPs.to_numpy().reshape(-1).astype(float)
pheno = pd.read_csv("step12a_pheno.txt", delimiter = "\t", header = None, usecols = [j])
pheno = pheno.to_numpy().reshape(-1)
E = pd.read_csv("step12a_Env_factor.txt", delimiter = "\t", header = None)
E = E.to_numpy().reshape(-1)
is_male = np.zeros(len(E))
pval = perm_test_bootstrap(SNP, pheno, E, is_male, 20000, 10000)
pd.DataFrame([pval]).to_csv("step0_GxE_p_values/" + str(j) + "_" + str(i) + ".txt", sep = "\t", header = False, index = False)


