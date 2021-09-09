import numpy as np
import pandas as pd
import argparse
import os
from tqdm import tqdm
from functools import reduce
from bed_reader import open_bed
from scipy.stats import linregress
from scipy.stats import yeojohnson as yj
import pdb

def nanlinregress(x, y):
    xval_indices = (np.isnan(x) == False)
    yval_indices = (np.isnan(y) == False)
    val_indices = np.logical_and(xval_indices, yval_indices)
    x, y = x[val_indices], y[val_indices]
    return(linregress(x, y))

parser = argparse.ArgumentParser()
parser.add_argument('--rsID', nargs = 1, type = str, action = "store", dest = "rsID")
args = parser.parse_args()
rsID = (args.rsID)[0]

hits = pd.read_csv("vQTL_hits.txt", delimiter = "\t")
hits["phenotype"] = hits.loc[:, "phenotype"].astype(int)
pheno, chr = hits.loc[hits["rsID"] == rsID, ["phenotype", "chr"]].values[0]

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
sorted_indices = sorted_main_indices[is_in_intersect]
new_fam = new_fam_main[sorted_indices] 
old_fam = pd.read_csv(old_fam_path, delim_whitespace = True, header = None)[0].to_numpy()
old_in_new = np.isin(old_fam, new_fam)
if np.all(old_fam == np.sort(old_fam)) == False:
    print("error1: code expects the fam file from step 4 to have sorted eids")
    exit()
if not np.all(new_fam == old_fam[old_in_new]):
    print("error2: code expects the fam file from step 4 to have sorted eids")
    exit()

geno_path = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr" + chr
pheno_path = "../step7_adjust_HF_for_covariates/phenotypes.txt"
rsIDs = pd.read_csv(geno_path  + ".bim", delim_whitespace = True, header = None)[1]
rsID_index = np.where(rsIDs == rsID)[0][0]
lb, ub = rsID_index - 2500, rsID_index + 2500
phenotypes = pd.read_csv(pheno_path, delimiter = "\t", header = None).loc[old_in_new, 1:].to_numpy()
phenotype = yj(phenotypes[:, pheno])[0]
genotypes = open_bed(geno_path  + ".bed", count_A1 = False, num_threads = 1).read(np.s_[sorted_indices, ])
missingness = np.sum(np.isnan(genotypes), axis = 0)/len(genotypes)
low_missingness = missingness < 0.02
keep_50_left = np.flip(np.cumsum(np.flip(low_missingness[:2500])))
keep_50_right = np.cumsum(low_missingness[2500:]) - 1
keep_100 = np.concatenate([keep_50_left, keep_50_right])
N = 1000
indices_to_keep = np.logical_and(low_missingness, keep_100 <= N)
genotypes = genotypes[:, indices_to_keep] 

main_effect = nanlinregress(genotypes[:, N], phenotype)[3]
if main_effect > 0.0001:
    print("main effect: " + str(main_effect))
    base_p_vals = []
    epistasis_p_vals1 = []
    epistasis_p_vals2 = []
    for i in range(len(genotypes[0])):
        base_p_vals.append(nanlinregress(genotypes[:, i], phenotype)[3])
        epistasis_p_vals1.append(nanlinregress(genotypes[:, N]*genotypes[:, i], phenotype)[3])
        epistasis_p_vals2.append(nanlinregress(genotypes[:, N]*(2 - genotypes[:, i]), phenotype)[3])

    print(np.min(epistasis_p_vals1))
    print(base_p_vals[np.argmin(epistasis_p_vals1)])
    print(np.min(epistasis_p_vals2))
    print(base_p_vals[np.argmin(epistasis_p_vals2)])
else:
    print("main effect: " + str(main_effect))