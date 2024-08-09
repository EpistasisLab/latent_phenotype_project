import numpy as np
import pandas as pd
import os
import pdb

'''
# formula for corrected m (g*) in https://pubmed.ncbi.nlm.nih.gov/22588983/
# formula for ICC: https://en.wikipedia.org/wiki/Intraclass_correlation
pheno_path = "../step9_regress_phenotypes_against_SNPs_PCA/QTL_phenotypes.txt"
X = pd.read_csv(pheno_path, delimiter = "\t", header = None).to_numpy()[:, :16]
N, K = len(X), len(X[0])
xn, x, s2 = np.mean(X, axis = 1), np.mean(X), np.var(X)
ICC = (K*np.sum((xn - x)**2)/(N*(K - 1)*s2)) - (1/(K - 1))
ICC = np.max([ICC, 0])
m = (K + 1) - (1 + (K - 1)*ICC)
'''

pb = 5E-8/16
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
main_effects = main_effects.drop_duplicates(["rsID", "pheno_index"])
Main = main_effects[["rsID", "p_null2", "pheno_index"]]
Main.to_csv("rsIDs_main_effects.txt", sep = "\t", header = True, index = False)
