import numpy as np
import pandas as pd
import argparse
import os
from functools import reduce
from bed_reader import open_bed
from scipy.stats import linregress
from scipy.stats import yeojohnson as yj
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
sorted_indices = sorted_main_indices[is_in_intersect]
new_fam = new_fam_main[sorted_indices] 
old_fam = pd.read_csv(old_fam_path, delim_whitespace = True, header = None)[0].to_numpy()
old_in_new = np.isin(old_fam, new_fam)
if not np.all(old_fam == np.sort(old_fam)):
    print("error1: code expects the fam file from step 4 to have sorted eids")
    exit()
if not np.all(new_fam == old_fam[old_in_new]):
    print("error2: code expects the fam file from step 4 to have sorted eids")
    exit()

geno_path = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr" + chr
pheno_path = "../step7_adjust_HF_for_covariates/phenotypes.txt"
rsIDs = pd.read_csv(geno_path  + ".bim", delim_whitespace = True, header = None)[1]
intervals = np.cumsum([0] + int(len(rsIDs)/100)*[100] + [len(rsIDs)%100])
phenotypes = pd.read_csv(pheno_path, delimiter = "\t", header = None).loc[old_in_new, 1:].to_numpy()
for i in range(16): phenotypes[:, i] = yj(phenotypes[:, i])[0]
if chr == '1' and pheno == 0:
    pd.DataFrame(phenotypes).to_csv("QTL_phenotypes.txt", sep = "\t", header = False, index = False)
phenotypes = phenotypes[:, pheno]
N = len(rsIDs)
out = np.zeros((N, 5))
i = 0
for j in tqdm(range(len(intervals) - 1)):
    lb, ub = intervals[j], intervals[j + 1]
    genotypes = open_bed(geno_path  + ".bed", count_A1 = False, num_threads = 1).read(np.s_[sorted_indices, lb:ub])
    missingness = np.sum(np.isnan(genotypes), axis = 0)/len(genotypes)
    valued_indices = np.isnan(genotypes) == False
    for k in range(len(genotypes[0])):
        g = genotypes[valued_indices[:, k], k]
        p = phenotypes[valued_indices[:, k]]
        if np.var(g) > 0 and np.var(p) > 0 and missingness[k] <= 0.02:
            out[i] += linregress(g, p)
            if out[i][3] < 5E-8/32:
                print(out[i][3])
        else:
            out[i] = np.array(5*[np.nan])
        i += 1

out = pd.DataFrame(out)
out.columns = ["slope", "intercept", "r", "p", "err"]
out["rsID"] = rsIDs
out = out[["rsID"] + ["slope", "intercept", "r", "p", "err"]]
out.to_csv("QTL_output/SNP_effects_chr" + chr + "_P" + str(pheno) + ".txt", sep = "\t", header = True, index = False)