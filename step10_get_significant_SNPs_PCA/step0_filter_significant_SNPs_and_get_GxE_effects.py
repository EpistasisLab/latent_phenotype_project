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
pheno_indices = list(range(16))

# maybe start with these: ['pack-years', 'annual-consumption', '874-average']
env_name = 'pack-years'
path = "../step7_adjust_HF_for_covariates_PCA/env_factors_for_step9.txt"
env_data = pd.read_csv(path, delimiter = "\t", usecols = ["eid", "22001-0.0", env_name])

QTL_pheno_path = "../step9_regress_phenotypes_against_SNPs_PCA/QTL_phenotypes.txt"
QTL_eids_path = "../step9_regress_phenotypes_against_SNPs_PCA/QTL_phenotypes_eids.txt"
QTL_phenotypes = pd.read_csv(QTL_pheno_path, delimiter = "\t", header = None).to_numpy()
QTL_eids = pd.read_csv(QTL_eids_path, delimiter = "\t", header = None)
QTL_eids.columns = ["eid"]

env_data = env_data.merge(QTL_eids, on = "eid", how = "inner")
is_male = env_data["22001-0.0"].to_numpy(dtype = "bool")
env_factor = env_data[env_name].to_numpy()

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
# new_fam_intersect = reduce(np.intersect1d, new_fams)
# all eids in QTL_eids exist in new_fam_intersec, but the converse is not true due to modifications. 
# the new new_fam_intersect term will serve the same function as previously. 
new_fam_intersect = QTL_eids["eid"].to_numpy()
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

def get_conditional_significance(genotypes_raw, phenotype_raw, E, rsIDs, top_rsID, best_hits, high_LD_hits, pb):

    genotypes = COPY(genotypes_raw).T
    for g in genotypes: g[np.isnan(g)] = np.nanmean(g)
    genotypes = genotypes.T
    
    phenotype = COPY(phenotype_raw)
    phenotype[np.isnan(phenotype)] = np.nanmean(phenotype)
    
    best_hits.append(top_rsID)
    top_indices = np.isin(rsIDs, best_hits)
    top_genotypes = genotypes[:, top_indices]
       
    g_top = genotypes[:, rsIDs == top_rsID].reshape(-1)
    corrs = np.array([nanlinregress(g_top, gi)[2] for gi in genotypes.T])
    high_LD_hits = np.concatenate([high_LD_hits, rsIDs[corrs**2 > 0.8]])
    exited_indices = np.logical_or(top_indices, np.isin(rsIDs, high_LD_hits))

    p_values = []
    for i in range(len(genotypes[0])):
        if exited_indices[i] == True:
            p_values.append(1)
        else:
            g = genotypes[:, i]
            G = np.concatenate([g.reshape(-1,1), top_genotypes], axis = 1)
            X = np.concatenate([G, E, np.ones((len(G), 1))], axis = 1)
            Z = COPY(X.T)
            for z in Z: z[np.isnan(z)] = np.nanmean(z)
            X_sub = np.concatenate([top_genotypes, E, np.ones((len(G), 1))], axis = 1)

            model = sm.OLS(phenotype, X)
            model_sub = sm.OLS(phenotype, X_sub)
            model_results = model.fit()
            model_sub_results = model_sub.fit()
            if model_results.llf == model_sub_results.llf:
                p_values.append(1)
            else:
                p_values.append(model_results.compare_lr_test(model_sub_results)[1])

    p_values = np.array(p_values)
    if np.any(p_values <= pb):
        
        kept_indices = np.logical_or((p_values <= pb), top_indices)
        next_genotypes = genotypes[:, kept_indices]
        next_rsIDs = np.array(rsIDs)[kept_indices]
        next_top_rsID = rsIDs[np.argmin(p_values)]
        next_args = [next_genotypes, phenotype, E, next_rsIDs, next_top_rsID, best_hits, high_LD_hits, pb]
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
        bim_file = pd.read_csv(path  + ".bim", delim_whitespace = True, header = None)
        rsIDs = bim_file[1]
        positions = bim_file[3]
        rsIDs_i_indices =  np.where(np.isin(rsIDs, QTLs_chr["rsID"]))[0]
        rsIDs_i = rsIDs[rsIDs_i_indices]
        positions_i = positions[rsIDs_i_indices].to_numpy()
        genotypes_getter = open_bed(path  + ".bed", count_A1 = False, num_threads = 1)
        genotypes = genotypes_getter.read(np.s_[sorted_indices,  rsIDs_i_indices])
        is_not_maf = np.nanmean(genotypes, axis = 0) > 1 
        genotypes[:, is_not_maf] = 2 - genotypes[:, is_not_maf]
        genotypes_enc = transform_genotypes(genotypes, QTL_P, env_factor)
        P = len(genotypes_enc[0])
        pb = 0.05/len(genotypes_enc[0])

        SL = 50
        segments = [SL*k for k in range(int(P/SL) + 1)]
        if P%SL != 0: segments += [SL*(int(P/SL)) + P%SL]
        best_rsIDs_i = []
        for k in tqdm(range(len(segments) - 1)):

            possible_hits_k = (rsIDs_i.to_numpy())[segments[k]:segments[k+1]]
            QTLs_chr_k = QTLs_chr.loc[QTLs_chr["rsID"].isin(possible_hits_k), :]
            if k == 0:
                prev_start = 0
                best_hits = []
                best_hits_innitial = COPY(best_hits)
                top_rsID = QTLs_chr_k["rsID"].to_numpy()[np.argmin(QTLs_chr_k["p_main"])]
            else:
                possible_start_positions = np.where(positions_i[segments[k] - 1] - positions_i > 1000000)[0]
                if len(possible_start_positions) == 0:
                    prev_start = 0
                else:
                    prev_start = np.max(possible_start_positions)
                best_hits = (rsIDs_i.to_numpy())[prev_start:segments[k]].tolist()
                best_hits_innitial = COPY(best_hits)
                if len(best_hits) == 0:
                    top_rsID = QTLs_chr_k["rsID"].to_numpy()[np.argmin(QTLs_chr_k["p_main"])]
                else:
                    top_rsID = best_hits.pop()

            rsIDs_ik = (rsIDs_i.to_numpy())[prev_start:segments[k+1]]
            genotypes_enc_k = genotypes_enc[:, prev_start:segments[k+1]]
            EF = env_factor.reshape(-1,1)
            best_rsIDs_ik = get_conditional_significance(genotypes_enc_k, QTL_P, EF, rsIDs_ik, top_rsID, best_hits, np.array([]), pb)
            best_rsIDs_i += np.setdiff1d(best_rsIDs_ik, best_hits_innitial).tolist()

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
