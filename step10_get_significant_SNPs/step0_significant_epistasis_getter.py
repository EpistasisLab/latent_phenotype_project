import numpy as np
import pandas as pd
import argparse
import os
from functools import reduce
from itertools import product
from bed_reader import open_bed
from scipy.stats import linregress
from scipy.stats import yeojohnson as yj
from scipy.stats import f_oneway
from copy import deepcopy as COPY
import pdb
from tqdm import tqdm

def nanlinregress(x, y):
    val_indices = np.logical_not(np.logical_or(np.isnan(x), np.isnan(y)))
    x, y = x[val_indices], y[val_indices]
    return(linregress(x, y))

if not os.path.exists("epistasis_output"):
    try: 
        os.mkdir("epistasis_output")
    except:
        pass

parser = argparse.ArgumentParser()
parser.add_argument('--chr', nargs = 1, type = str, action = "store", dest = "chr")
parser.add_argument('--rsID', nargs = 1, type = str, action = "store", dest = "rsID")
parser.add_argument('--pheno', nargs = 1, type = str, action = "store", dest = "pheno")
args = parser.parse_args()
rsID = (args.rsID)[0]
chr = (args.chr)[0]
phenotype_index = (args.pheno)[0]

vQTL_info = pd.read_csv("vQTL_hits_female.txt", delimiter = "\t", dtype = str)
cond = np.logical_and(vQTL_info["rsID"] == rsID, vQTL_info["phenotype"] == phenotype_index)
alpha = vQTL_info.loc[cond, "alpha"].item()
rsID_info = vQTL_info[cond]
chromosomes = np.arange(1, 23).astype(str).tolist() + ["X", "XY", "Y", "MT"]
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
old_fam = pd.read_csv(old_fam_path, delim_whitespace = True, header = None)[0].to_numpy()
old_in_new = np.isin(old_fam, new_fam)
if np.all(old_fam == np.sort(old_fam)) == False:
    print("error2: code expects the fam file from step 4 to have sorted eids")
    exit()
if not np.all(new_fam == old_fam[old_in_new]):
    print("error3: code expects the fam file from step 4 to have sorted eids")
    exit()

# allows for seperate regression of males and females
path = "../step7_adjust_HF_for_covariates/X.txt"
gender_data = pd.read_csv(path, delimiter = "\t", usecols = ["eid", "22001-0.0"], dtype = int)
new_fam_df = pd.DataFrame(np.array([new_fam, np.arange(len(new_fam))]).T)
new_fam_df.columns = ["eid", "index"] 
gender_data = gender_data.merge(new_fam_df, on = "eid", how = "inner").sort_values(by = "index")
is_male = gender_data[ "22001-0.0"].to_numpy(dtype = "bool")

geno_path = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr" + chr
pheno_path = "../step7_adjust_HF_for_covariates/phenotypes_cleaned_data.txt"
phenotypes = pd.read_csv(pheno_path, delimiter = "\t")
good_indices = np.isin(phenotypes["eid"].to_numpy(), new_fam)
phenotype = phenotypes.loc[good_indices, phenotype_index[0]]
phenotype = yj(phenotype.to_numpy())[0]

vQTL_chromosome = rsID_info["chr"].to_numpy()[0]
vQTL_path = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr" + str(vQTL_chromosome) 
new_fam_vQTL = pd.read_csv(vQTL_path + ".fam", delim_whitespace = True, header = None)[1].to_numpy()
is_in_intersect_VQTL = np.isin(new_fam_vQTL, new_fam_intersect)
sorted_main_indices_vQTL = np.argsort(new_fam_vQTL)
sorted_indices_vQTL = sorted_main_indices_vQTL[is_in_intersect_VQTL[sorted_main_indices_vQTL]]
vQTL_chromosome_rsIDs = pd.read_csv(vQTL_path  + ".bim", delim_whitespace = True, header = None) 
vQTL_pos = np.where(vQTL_chromosome_rsIDs == rsID)[0][0]
filter = np.s_[sorted_indices_vQTL, vQTL_pos]
vQTL_genotypes = open_bed(vQTL_path  + ".bed", count_A1 = False, num_threads = 1).read(filter).reshape(-1)
cond_partial = np.logical_and(vQTL_info["rsID"] == rsID, vQTL_info["chr"] == chr)
cond = np.logical_and(cond_partial, vQTL_info["phenotype"] == phenotype_index)
alpha1 = rsID_info["alpha"].to_numpy(dtype = float)[0]
vQTL_genotypes[vQTL_genotypes == 1] = alpha1
vQTL_genotypes[vQTL_genotypes == 2] = 1

# SNPs examined for QTLs and VQTLs are the same.
prefix = "../step9_regress_phenotypes_against_SNPs/QTL_output/"
partial_path = prefix + "QTL_effects_chr" + chr + "_P" + phenotype_index[0]
outs_gender = []
prefix = "../step9_regress_phenotypes_against_SNPs/QTL_output/"
rsIDs = []
for gender, gender_indices in zip(["male", "female"], [is_male, is_male == False]):
    if chr == "Y" and gender == "female":
        continue

    phenotype_gender = phenotype[gender_indices]
    vQTL_genotypes_gender = vQTL_genotypes[gender_indices]
    m, b, void, void, void = nanlinregress(vQTL_genotypes_gender, phenotype_gender)
    residuals1 = phenotype_gender - (m*vQTL_genotypes_gender + b)

    path = partial_path + "_" + gender + ".txt"
    gender_geno_info = pd.read_csv(path, delimiter = "\t")
    geno_bim_path = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr" + chr + ".bim"
    geno_bim_rsIDs = pd.read_csv(geno_bim_path, delim_whitespace = True, header = None, usecols = [1])
    # good meaning population missingness =< 0.02 and maf >= 0.005, as previously filtered
    is_good_RSID = geno_bim_rsIDs.isin(gender_geno_info.dropna()["rsID"]).to_numpy().reshape(-1)
    rsIDs.append(geno_bim_rsIDs[is_good_RSID].to_numpy().reshape(-1))
    # checks that the different bim sources are ordered in the same way
    source1 = gender_geno_info.dropna()["rsID"].to_numpy().reshape(-1)
    source2 = geno_bim_rsIDs[is_good_RSID].to_numpy().reshape(-1)
    if np.any(source1 != source2):
        print("exiting: the bim file's rsIDs are not ordered the same as the step8 analysis file's rsIDs")
        exit()
    info_packets = gender_geno_info.dropna()
    good_rsID_indices = np.where(is_good_RSID)[0]
    N = np.sum(is_good_RSID)
    intervals = np.cumsum([0] + int(N/100)*[100] + [N%100])
    out = np.zeros(N)
    i = 0
    for j in tqdm(range(len(intervals) - 1)):
        lb, ub = intervals[j], intervals[j + 1]
        interval_indices = good_rsID_indices[lb:ub]
        filter = np.s_[sorted_indices[gender_indices], interval_indices]
        genotypes = open_bed(geno_path  + ".bed", count_A1 = False, num_threads = 1).read(filter)
        info_subpackets = info_packets.loc[interval_indices, :]
        for k in range(len(genotypes[0])):
            alpha2 = info_subpackets[["alpha"]].to_numpy()[k][0]
            g = COPY(genotypes[:, k])
            if alpha2 == -1:
                vals, counts = np.unique(g[np.isnan(g) == False], return_counts = True)
                if len(vals) == 3:
                    val_to_change = vals[np.argmin(counts)]
                    g[g ==  val_to_change] = 1
            else:
                g[g == 1] = alpha2
                g[g == 2] = 1
            m2, b2, v, v, v = nanlinregress(g, residuals1)
            residuals2 = residuals1 - (m2*g + b2)
            data = np.array([vQTL_genotypes_gender, g, residuals2]).T
            data = data[np.any(np.isnan(data), axis = 1) == False]
            X, Y = np.round(data[:, 0:2], 6), data[:, 2]
            X_hashes = pd.util.hash_pandas_object(pd.DataFrame(X), index = False).to_numpy()
            hash_IDs = np.unique(X_hashes)
            categories = X_hashes == hash_IDs.reshape(-1, 1)
            grouped_phenotypes = []
            for ind in categories:
                phenotype_i = Y[ind]
                if len(phenotype_i) > 200:
                    if len(phenotype_i) < 1000:
                        p_low = phenotype_i > np.percentile(phenotype_i, 2.5)
                        p_high = phenotype_i < np.percentile(phenotype_i, 97.5)
                        phenotype_i = phenotype_i[np.logical_and(p_low, p_high)]
                    grouped_phenotypes.append(phenotype_i)
            out[i] += f_oneway(*grouped_phenotypes)[1]
            i += 1
    outs_gender.append(out)

pdb.set_trace()
for i, out, gender in zip([0,1], outs_gender, ["_male", "_female"]):
    out = pd.DataFrame(out)
    out.columns = ["slope", "intercept", "r", "p", "err"]
    out["rsID"] = rsIDs[i]
    out = out[["rsID"] + ["slope", "intercept", "r", "p", "err"]]
    path = "epistasis_output/SNP_effects_chr" + chr + "_" + rsID + "_P" + phenotype_index[0] + gender + ".txt"
    out.to_csv(path, sep = "\t", header = True, index = False)