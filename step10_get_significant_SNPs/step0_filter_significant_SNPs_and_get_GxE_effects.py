import numpy as np
import pandas as pd
import argparse
import os
from functools import reduce
from tqdm import tqdm
from copy import deepcopy as COPY
from scipy.stats import linregress
import statsmodels.api as sm
from scipy.stats import chi2
from bed_reader import open_bed
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
import pdb
import matplotlib as mpl
from scipy.stats import norm
mpl.rcParams['agg.path.chunksize'] = 10000
from statsmodels.stats.outliers_influence import summary_table

parser = argparse.ArgumentParser()
parser.add_argument('--chr', nargs = 1, type = int, action = "store", dest = "chr")
parser.add_argument('--name', nargs = 1, type = str, action = "store", dest = "name")
args = parser.parse_args()
chr_index = (args.chr)[0] - 1
name = (args.name)[0]

if not os.path.exists("hits_QTL" + name):
    os.mkdir("hits_QTL" + name)
if not os.path.exists("hits_GxE_p_vals_getters" + name):
    os.mkdir("hits_GxE_p_vals_getters" + name)

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]

all_chromosomes = COPY(chromosomes)
pheno_indices = list(range(17))

# maybe start with these: ['pack-years', 'annual-consumption', '874-average']
env_name = 'annual-consumption'
path = "../step9_regress_phenotypes_against_SNPs/env_factors_cleaned.txt"
env_data = pd.read_csv(path, delimiter = "\t", usecols = ["22001-0.0", env_name])
is_male = env_data["22001-0.0"].to_numpy(dtype = "bool")
env_factor = env_data[env_name].to_numpy()
env_factor[np.isnan(env_factor)] = np.nanmedian(env_factor)

QTL_pheno_path = "../step9_regress_phenotypes_against_SNPs/QTL_phenotypes.txt"
QTL_phenotypes = pd.read_csv(QTL_pheno_path, delimiter = "\t", header = None).to_numpy()

QTL_paths = ["significant_QTLs" + name + "/significant_QTLs_P" + str(i) + ".txt" for i in pheno_indices]
QTLs_info = [pd.read_csv(path, delimiter = "\t") for path in QTL_paths]
for info in QTLs_info: info["chr"] = info["chr"].astype(str)

chr_plink_paths = ["significant_SNPs_plink_files" + name + "/significant_SNPs_chr" + chr for chr in chromosomes]  
chr_has_hits = []
for path in chr_plink_paths:
    if os.path.isfile(path + ".bed"):
        chr_has_hits.append(True)
    else:
        chr_has_hits.append(False)
chr_has_hits = np.array(chr_has_hits)
chromosomes = np.array(chromosomes)[chr_has_hits]
chr_plink_paths = np.array(chr_plink_paths)[chr_has_hits]

new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
new_fam_paths = np.array([new_fam_path_prefix + i + ".fam" for i in all_chromosomes])
old_fam_path = "../step4_remove_relatives/UKB_samples_unrelated.fam"
new_fams = [pd.read_csv(path, delim_whitespace = True, header = None)[1].to_numpy() for path in new_fam_paths]
new_fam_intersect = reduce(np.intersect1d, new_fams)
new_fam_paths = new_fam_paths[chr_has_hits]

def nanlinregress(x, y):
    xval_indices = (np.isnan(x) == False)
    yval_indices = (np.isnan(y) == False)
    val_indices = np.logical_and(xval_indices, yval_indices)
    x, y = x[val_indices], y[val_indices]
    return(linregress(x, y))


def transform_genotypes(genotypes, phenotypes, Env_factor):
    
    for i in range(len(genotypes[0])):

        g0 = genotypes[:, i]
        valued_indices = np.logical_or(np.isnan(phenotypes), np.isnan(g0)) == False
        g = genotypes[valued_indices, i]
        p = phenotypes[valued_indices]
        E1 = Env_factor[valued_indices]

        g2_vals = (g == 2)
        if np.all(np.sum(g2_vals, axis = 0) > 1000) and chr != "Y": 

            g1_vals = (g == 1)
            X = np.array([g1_vals, g2_vals, (g1_vals)*E1, (g2_vals)*E1, E1]).T
            X2 = np.array([g1_vals, g2_vals]).T
            y = (p - np.mean(p))

            Betas = LinearRegression().fit(X, y).coef_
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

            Betas = LinearRegression().fit(X, y).coef_
            enc = Betas[0] + Betas[1]*E1[g1_vals]  
            g_old = COPY(g)
            g[g1_vals] = enc 

        genotypes[valued_indices, i] = g

    return(genotypes)

def get_conditional_significance(genotypes, phenotype, E, rsIDs, top_rsID, best_hits):

    k = len(top_rsID)
    pb = 0.05/(len(rsIDs))
    best_hits.append(top_rsID)
    top_indices = np.isin(rsIDs, best_hits)
    top_genotypes = genotypes[:, top_indices]
    p_values = []
    for i in range(len(genotypes[0])):
        if top_indices[i] == True:
            p_values.append(1)
        else:
            g = genotypes[:, i]
            G = np.concatenate([g.reshape(-1,1), top_genotypes], axis = 1)
            X = np.concatenate([G, E, np.ones((len(G), 1))], axis = 1)
            X += np.random.normal(0, 1E-6, X.shape)
            X_sub = np.concatenate([top_genotypes, E, np.ones((len(G), 1))], axis = 1)
            X_sub += np.random.normal(0, 1E-6, X_sub.shape)
            vals1 = np.isnan(phenotype) == False
            vals2 = np.all(np.isnan(G) == False, axis = 1)
            vals3 = np.all(np.isnan(top_genotypes) == False, axis = 1)
            vals4 = np.logical_and(vals1, vals2)
            vals5 = np.logical_and(vals1, vals3)
            model = sm.OLS(phenotype[vals4], X[vals4])
            model_sub = sm.OLS(phenotype[vals5], X_sub[vals5])
            model_results = model.fit()
            model_sub_results = model_sub.fit()
            p_values.append(model_results.compare_lr_test(model_sub_results)[1])
    p_values = np.array(p_values)
    if np.any(p_values <= pb):
        kept_indices = np.logical_or((p_values <= pb), top_indices)
        next_genotypes = genotypes[:, kept_indices]
        next_rsIDs = np.array(rsIDs)[kept_indices]
        next_top_rsID = rsIDs[np.argmin(p_values)]
        next_args = [next_genotypes, phenotype, E, next_rsIDs, next_top_rsID, best_hits]
        return(get_conditional_significance(*next_args))
    else:
        return(best_hits)

best_rsIDs = []
QTLs_chr_all = []
for i in pheno_indices:
    if i == 16 and chromosomes[chr_index] == "6":
        continue

    print("i: " + str(i))
    QTL_P = QTL_phenotypes[:, i]
    QTLs = QTLs_info[i]

    chr = chromosomes[chr_index]
    path = chr_plink_paths[chr_index]
    print(chr)       
 
    main_path = new_fam_paths[chr_index]
    new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
    is_in_intersect = np.isin(new_fam_main, new_fam_intersect)
    sorted_main_indices = np.argsort(new_fam_main)
    sorted_indices = sorted_main_indices[is_in_intersect[sorted_main_indices]]

    QTLs_chr = QTLs.loc[QTLs["chr"] == chr, :]
    QTLs_chr["pheno_index"] = i
    QTLs_chr_all.append(QTLs_chr)
    if len(QTLs_chr) > 0:
        top_rsID = QTLs_chr["rsID"].to_numpy()[np.argmin(QTLs_chr["p_main"])]
        rsIDs = pd.read_csv(path  + ".bim", delim_whitespace = True, header = None)[1]
        rsIDs_i_indices =  np.where(np.isin(rsIDs, QTLs_chr["rsID"]))[0]
        rsIDs_i = rsIDs[rsIDs_i_indices]
        genotypes_getter = open_bed(path  + ".bed", count_A1 = False, num_threads = 1)
        genotypes = genotypes_getter.read(np.s_[sorted_indices,  rsIDs_i_indices])
        is_not_maf = np.nanmean(genotypes, axis = 0) > 1 
        genotypes[:, is_not_maf] = 2 - genotypes[:, is_not_maf]
        genotypes_enc = transform_genotypes(genotypes, QTL_P, env_factor)
        
        best_rsIDs_i = get_conditional_significance(genotypes_enc, QTL_P, env_factor.reshape(-1,1), rsIDs_i.to_numpy(), top_rsID, [])
        best_rsIDs_i_loc = np.isin(rsIDs_i.to_numpy(), best_rsIDs_i)
        L = len(best_rsIDs_i)
        best_rsID_geno = genotypes[:, best_rsIDs_i_loc].T
        for k in range(L): best_rsID_geno[k, np.isnan(best_rsID_geno[k, :])] = np.nanmean(best_rsID_geno[k, :])
        max_R2 = np.max((np.corrcoef(best_rsID_geno)**2 - np.eye(L)))
        print(max_R2)
        best_rsIDs.append(np.array([best_rsIDs_i, L*[i], L*[max_R2]]).T)

QTLs_chr_all = pd.concat(QTLs_chr_all)
best_rsIDs2 = pd.DataFrame(np.concatenate(best_rsIDs, axis = 0))
best_rsIDs2.columns = ["rsID", "pheno_index", "max_LD"]
best_rsIDs2["pheno_index"] = best_rsIDs2["pheno_index"].astype(int)
best_rsIDs3 = best_rsIDs2.merge(QTLs_chr_all, on = ["rsID", "pheno_index"], how = "inner")
path = "hits_QTL" + name + "/QTL_hits_chr" + chromosomes[chr_index] + ".txt"
best_rsIDs3.to_csv(path, sep = "\t", header = True, index = False)
best_rsIDs3.index = np.arange(len(best_rsIDs3))

for i in range(len(best_rsIDs3)):
    info_i = best_rsIDs3.loc[i, ["rsID", "pheno_index", "chr"]].to_numpy(dtype = str)
    name_parts = COPY(info_i)
    name_parts[0] =  name_parts[0].split(":")[-1]
    file = open("hits_GxE_p_vals_getters" + name + "/" + "_".join(name_parts) + ".sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#BSUB -J GxE_p_getter\n")
    file.write("#BSUB -o hits_GxE_p_vals_getters" + name + "/GxE_p_getters.out\n")
    file.write("#BSUB -e hits_GxE_p_vals_getters" + name + "/GxE_p_getters.err\n")
    file.write("source activate torch_env2\n\n")
    file.write("python step0_compute_GxE_p_values.py ")
    file.write("--rsID " + info_i[0] + " --pheno_index " + str(info_i[1]) + " --chr " + info_i[2] + " --name " + name)
    file.close()


print(1)
