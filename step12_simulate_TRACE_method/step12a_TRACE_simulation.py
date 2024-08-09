import numpy as np
import pandas as pd
import os
from scipy.stats import pearsonr
from scipy.stats import logistic
from scipy.stats import uniform
from scipy.stats import kstest
from sklearn.decomposition import PCA
import statsmodels.api as sm
from matplotlib import pyplot as plt
import pdb

def standardize(x): return((x - np.mean(x))/np.std(x))

N, P, L, Z, J = 200000, 100, 14, 4, 1800
if not os.path.exists("step12a_SNPs.txt"):

    # simulating causal SNPs and their minor allele frequencies (mafs)
    main_eff_mafs = np.random.uniform(0.05, 0.5, P)
    GxE_eff_mafs = np.random.uniform(0.05, 0.5, P)
    main_eff_SNPs = np.random.binomial(n = 2, p = main_eff_mafs, size = (N, P))
    GxE_eff_SNPs = np.random.binomial(n = 2, p = GxE_eff_mafs, size = (N, P))

    # simulating logistic liability components from SNP associations
    SNP_assign_ubs = np.round(np.cumsum(np.repeat(P/L, L)),0).astype(int)
    SNP_assign_lbs = np.concatenate([np.zeros(1), SNP_assign_ubs[:-1]]).astype(int)
    liability_components = np.zeros((L, N))
    env_factor = np.random.normal(0, 1, N)
    # assigns main eff, GxE eff, E eff, and random noise
    # TODO: make env effect comperable to genetic effect

    h_narrow, h_GxE  = [], []
    # true heritability will be 20% main + 40% GxE + 40% noise
    # true GxE is statistically E in part. 
    for i, lb, ub in zip(np.arange(L), SNP_assign_lbs, SNP_assign_ubs):
        sig_main, sig_GxE, sig_noise = (8**(1/2)), (16**(1/2)), (76**(1/2))
        main_effs = main_eff_SNPs[:, lb:ub]
        main_component = np.sum(main_effs, axis = 1)
        main_component = standardize(main_component)*sig_main
        GxE_effs = GxE_eff_SNPs[:, lb:ub]*env_factor.reshape(-1,1)
        GxE_component = np.sum(GxE_effs, axis = 1)
        GxE_component = standardize(GxE_component)*sig_GxE
        noise_component = np.random.normal(0, sig_noise, N)
        liability_components[i] += main_component
        liability_components[i] += GxE_component
        liability_components[i] += noise_component        
    
        '''
        # This computes the heritability of the underlying continuous phenotypes
        X = np.concatenate([main_effs, np.ones((N, 1))], axis = 1)
        model = sm.OLS(liability_components[i], X).fit()
        h_narrow.append(model.rsquared)
    
        X = np.concatenate([main_effs, GxE_effs, np.ones((N, 1))], axis = 1)
        model = sm.OLS(liability_components[i], X).fit()
        h_GxE.append(model.rsquared)
        '''

    weights = [[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]]*20 
    weights += [[0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]]*20 
    weights += [[0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]]*20 
    weights += [[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0]]*20
    weights += [[0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]]*20
    weights = np.array(weights)

    sums = np.matmul(weights, liability_components)
    exponent = (sums - np.mean(sums))/np.std(sums) - 2
    probs = logistic.cdf(exponent)
    ICD_codes = np.random.binomial(n = 1, p = probs, size = probs.shape)

    junk_mafs = np.random.uniform(0.05, 0.5, J)
    junk_SNPs = np.random.binomial(n = 2, p = junk_mafs, size = (N, J))
    X = np.concatenate([main_eff_SNPs, GxE_eff_SNPs, junk_SNPs], axis = 1)
 
    pca = PCA(n_components = Z)
    pca.fit(ICD_codes.T)
    PCs = pca.transform(ICD_codes.T)
    loadings = pca.components_
    ICD_reconstructed = np.matmul(PCs, loadings) + np.mean(ICD_codes.T, axis = 0)
    residuals = ICD_codes[0] - ICD_reconstructed[:, 0]
    Y = np.concatenate([PCs, residuals.reshape(-1,1)], axis = 1)

    pdb.set_trace()
    pd.DataFrame(ICD_codes.T).to_csv("step12a_ICD_codes.txt", sep = "\t", header = False, index = False)
    pd.DataFrame(Y).to_csv("step12a_pheno.txt", sep = "\t", header = False, index = False)
    pd.DataFrame(X.astype(int)).to_csv("step12a_SNPs.txt", sep = "\t", header = False, index = False)
    pd.DataFrame(env_factor.reshape(-1,1)).to_csv("step12a_Env_factor.txt", sep = "\t", header = False, index = False)
else:
    X = pd.read_csv("step12a_SNPs.txt", delimiter = "\t", header = None).to_numpy(dtype = float)
    Y = pd.read_csv("step12a_pheno.txt", delimiter = "\t", header = None).to_numpy()
    ICD_codes = pd.read_csv("step12a_ICD_codes.txt", delimiter = "\t", header = None).to_numpy().T
    env_factor = pd.read_csv("step12a_Env_factor.txt", delimiter = "\t", header = None)[0].to_numpy()

E = env_factor.reshape(-1,1)
GxE = X*E
con = np.ones((len(E), 1))
data = np.concatenate([X, GxE, E, con], axis = 1)
Q = 2*P+J

normal_main_eff_pvals = []
normal_main_eff_null_pvals = []
normal_GxE_eff_pvals = []
normal_GxE_eff_null_pvals = []
for k in [0, 20, 40, 60, 80]:
    model = sm.Logit(ICD_codes[k], data).fit()
    normal_main_eff_pvals.append(model.pvalues[0:P])
    normal_main_eff_null_pvals.append(model.pvalues[P:Q])
    normal_GxE_eff_pvals.append(model.pvalues[(Q + P):(Q + 2*P)])
    Sa, Sb = model.pvalues[Q:(Q + P)], model.pvalues[(Q + 2*P):2*Q]
    next_set = np.concatenate([Sa, Sb])
    normal_GxE_eff_null_pvals.append(next_set)

df = pd.DataFrame(np.array(normal_main_eff_pvals).T)
df.to_csv("step12a_normal_main_eff_pvals.txt", sep = "\t", header = False, index = False)
df = pd.DataFrame(np.array(normal_main_eff_null_pvals).T)
df.to_csv("step12a_normal_main_eff_null_pvals.txt", sep = "\t", header = False, index = False)
df = pd.DataFrame(np.array(normal_GxE_eff_pvals).T)
df.to_csv("step12a_normal_GxE_eff_pvals.txt", sep = "\t", header = False, index = False)
df = pd.DataFrame(np.array(normal_GxE_eff_null_pvals).T)
df.to_csv("step12a_normal_GxE_eff_null_pvals.txt", sep = "\t", header = False, index = False)

TRACE_main_eff_pvals = []
TRACE_main_eff_null_pvals = []
data = np.concatenate([X, con], axis = 1)
for i in range(Z + 1):

    model = sm.OLS(Y[:, i], data).fit()
    TRACE_main_eff_pvals.append(model.pvalues[0:P])
    TRACE_main_eff_null_pvals.append(model.pvalues[2*P:-1])

pd.DataFrame(TRACE_main_eff_pvals).to_csv("step12a_TRACE_main_eff_pvals.txt", sep = "\t", header = False, index = False)
pd.DataFrame(TRACE_main_eff_null_pvals).to_csv("step12a_TRACE_main_eff_null_pvals.txt", sep = "\t", header = False, index = False)