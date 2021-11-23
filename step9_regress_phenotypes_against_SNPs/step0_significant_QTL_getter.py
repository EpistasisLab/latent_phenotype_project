import numpy as np
import pandas as pd
import argparse
import os
from functools import reduce
from copy import deepcopy as COPY
from bed_reader import open_bed
from scipy.stats import linregress
from scipy.stats import mode
from scipy.stats import yeojohnson as yj
from sklearn.linear_model import LinearRegression
import pdb
from tqdm import tqdm

if not os.path.exists("QTL_output"):
    try: 
        os.mkdir("QTL_output")
    except:
        pass

parser = argparse.ArgumentParser()
parser.add_argument('--chr', nargs = 1, type = str, action = "store", dest = "chr")
parser.add_argument('--pheno', nargs = 1, type = int, action = "store", dest = "pheno")
args = parser.parse_args()
chr = (args.chr)[0]
pheno = (args.pheno)[0]

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
new_fam_paths = [new_fam_path_prefix + i + ".fam" for i in chromosomes]
old_fam_path = "../step4_remove_relatives/UKB_samples_unrelated.fam"
new_fams = [pd.read_csv(path, delim_whitespace = True, header = None)[1].to_numpy() for path in new_fam_paths]
main_path = np.array(new_fam_paths)[np.array(chromosomes) == chr][0]
new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
new_fam_intersect = reduce(np.intersect1d, new_fams)
is_in_intersect = np.isin(new_fam_main, new_fam_intersect)
sorted_main_indices = np.argsort(new_fam_main)
sorted_indices = sorted_main_indices[is_in_intersect[sorted_main_indices]]
new_fam = new_fam_main[sorted_indices]

# allows for seperate regression of males and females
path = "../step7_adjust_HF_for_covariates/X.txt"
gender_data = pd.read_csv(path, delimiter = "\t", usecols = ["eid", "22001-0.0"], dtype = int)
new_fam_df = pd.DataFrame(np.array([new_fam, np.arange(len(new_fam))]).T)
new_fam_df.columns = ["eid", "index"] 
gender_data = gender_data.merge(new_fam_df, on = "eid", how = "inner").sort_values(by = "index")
is_male = gender_data[ "22001-0.0"].to_numpy(dtype = "bool")

old_fam = pd.read_csv(old_fam_path, delim_whitespace = True, header = None)[0].to_numpy()
old_in_new = np.isin(old_fam, new_fam)
if not np.all(old_fam == np.sort(old_fam)):
    print("error1: code expects the fam file from step 4 to have sorted eids")
    exit()
if not np.all(new_fam == old_fam[old_in_new]):
    print("error2: code expects the fam file from step 4 to have sorted eids")
    exit()

geno_path = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr" + chr
pheno_path = "../step7_adjust_HF_for_covariates/phenotypes_cleaned_data.txt"
rsIDs = pd.read_csv(geno_path  + ".bim", delim_whitespace = True, header = None)[1]
intervals = np.cumsum([0] + int(len(rsIDs)/100)*[100] + [len(rsIDs)%100])
phenotypes = pd.read_csv(pheno_path, delimiter = "\t").loc[old_in_new, np.arange(16).astype(str)].to_numpy()
for i in range(16): 
    val_indices = np.isnan(phenotypes[:, i]) == False
    phenotypes[val_indices, i] = yj(phenotypes[val_indices, i])[0]
if chr == '1' and pheno == 0:
    pd.DataFrame(phenotypes).to_csv("QTL_phenotypes.txt", sep = "\t", header = False, index = False)
    pd.DataFrame(is_male).to_csv("is_male.txt", sep = "\t", header = False, index = False)
    env_factors_path = "../step7_adjust_HF_for_covariates/env_factors.txt"
    env_factors = pd.read_csv(env_factors_path, delimiter = "\t", header = None).loc[old_in_new, 1:]
    env_factors.to_csv("env_factors.txt", sep = "\t", header = False, index = False)
phenotypes = phenotypes[:, pheno]
N = len(rsIDs)
out = [np.zeros((N, 7)), np.zeros((N, 7))]
i = 0
for j in tqdm(range(len(intervals) - 1)):

    lb, ub = intervals[j], intervals[j + 1]
    genotypes = open_bed(geno_path  + ".bed", count_A1 = False, num_threads = 1).read(np.s_[sorted_indices, lb:ub])
    if chr == "Y":
        missingness = np.sum(np.isnan(genotypes[is_male]), axis = 0)/len(genotypes[is_male])
    else:
        missingness = np.sum(np.isnan(genotypes), axis = 0)/len(genotypes)
    afs = np.nansum(genotypes, axis = 0)/(np.sum(np.isnan(genotypes) == False, axis = 0)*2)
    af_is_not_maf = afs > 1 - afs
    genotypes[:, af_is_not_maf] = 2 - genotypes[:, af_is_not_maf]
    valued_indices = np.isnan(genotypes) == False

    for k in range(len(genotypes[0])):
        g_all = genotypes[valued_indices[:, k], k]
        p_all = phenotypes[valued_indices[:, k]]
        gender = is_male[valued_indices[:, k]]
        g_gender = [g_all[gender == 1], g_all[gender == 0]]
        p_gender = [p_all[gender == 1], p_all[gender == 0]]
        for out_index, g, p in zip([0,1], g_gender, p_gender):
            valued_indices2 = np.isnan(p) == False
            g = g[valued_indices2]
            p = p[valued_indices2]
            if chr == "Y" and out_index == 1:
                continue
            maf = np.sum(g)/(2*len(g))
            if np.var(p) > 0 and missingness[k] <= 0.02 and maf >= 0.005:
                X = np.array([g == 1, g == 2]).T
                y = (p - np.mean(p))
                if np.all(np.sum(X, axis = 0) > 1000): 
                    Betas = LinearRegression().fit(X, y).coef_
                    out[out_index][i, 6] = Betas[0]/Betas[1]
                    g[g == 1] = Betas[0]/Betas[1]
                    g[g == 2] = 1                
                else:
                    out[out_index][i, 6] = -1
                    g[g == 2] = 1
                out[out_index][i, 1:6] = linregress(g, p)
            else:
                out[out_index][i] = np.array(7*[np.nan])
        i += 1

for o, label in zip(out, ["_male", "_female"]):
    if chr == "Y" and label == "_female":
        continue
    o_df = pd.DataFrame(o)
    o_df.columns = ["rsID", "slope", "intercept", "r", "p", "err", "alpha"]
    o_df["rsID"] = rsIDs
    o_df.to_csv("QTL_output/QTL_effects_chr" + chr + "_P" + str(pheno) + label + ".txt", sep = "\t", header = True, index = False)