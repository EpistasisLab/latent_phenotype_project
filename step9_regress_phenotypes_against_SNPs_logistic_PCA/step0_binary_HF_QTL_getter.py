import numpy as np
import pandas as pd
import argparse
import os
from functools import reduce
from copy import deepcopy as COPY
from bed_reader import open_bed
import statsmodels.api as sm
from scipy.stats import chi2
import pdb
from tqdm import tqdm

if not os.path.exists("binary_HF_QTL_output"):
    try: 
        os.mkdir("binary_HF_QTL_output")
    except:
        pass

parser = argparse.ArgumentParser()
parser.add_argument('--chr', nargs = 1, type = str, action = "store", dest = "chr")
parser.add_argument('--pheno', nargs = 1, type = str, action = "store", dest = "pheno")
parser.add_argument('--jstart', nargs = 1, type = int, action = "store", dest = "jstart")
args = parser.parse_args()
chr = (args.chr)[0]
pheno = (args.pheno)[0]
j_start = (args.jstart)[0]

suffix = "_chr" + chr + "_P" + str(pheno) + "_jstart" + str(j_start)
filename = "binary_HF_QTL_output/binary_HF_QTL_effects" + suffix + ".txt"
if os.path.exists(filename):
    exit()

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
new_fam_paths = [new_fam_path_prefix + i + ".fam" for i in chromosomes]
new_fams = [pd.read_csv(path, delim_whitespace = True, header = None)[1].to_numpy() for path in new_fam_paths]
main_path = np.array(new_fam_paths)[np.array(chromosomes) == chr][0]
new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
new_fam_intersect = reduce(np.intersect1d, new_fams)
is_in_intersect = np.isin(new_fam_main, new_fam_intersect)
sorted_main_indices = np.argsort(new_fam_main)
sorted_indices = sorted_main_indices[is_in_intersect[sorted_main_indices]]
new_fam = new_fam_main[sorted_indices]
new_fam_df = pd.DataFrame(np.array([new_fam, np.arange(len(new_fam))]).T)
new_fam_df.columns = ["eid", "index"] 

old_fam_path = "../step4_remove_relatives/UKB_samples_unrelated.fam"
old_fam = pd.read_csv(old_fam_path, delim_whitespace = True, header = None)[0].to_numpy()
old_in_new = np.isin(old_fam, new_fam)
# changing the definition of heart failure led to a slight change in the eids, which cannot be reflected in step 8
# only unrelated individuals under the previous procedure (step 8) and the current procedure (step 4) are included
# This does include all all-cause heart failure cases and only removes a few thousand controls.
new_in_old = np.isin(new_fam, old_fam)
if not np.all(old_fam == np.sort(old_fam)):
    print("error1: code expects the fam file from step 4 to have sorted eids")
    exit()
if not np.all(new_fam[new_in_old] == old_fam[old_in_new]):
    print("error2: code expects the fam file from step 4 to have sorted eids")
    exit()

geno_path = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr" + chr
pheno_path = "../step7_adjust_HF_for_covariates_logistic_PCA/binary_HF_values.txt"
PCs_path = "../step7_adjust_HF_for_covariates_logistic_PCA/binary_HF_genetic_PCs.txt"
Betas_path = "../step7_adjust_HF_for_covariates_logistic_PCA/binary_HF_PC_betas.txt"
rsIDs = pd.read_csv(geno_path  + ".bim", delim_whitespace = True, header = None)[1]

phenotypes_info = pd.read_csv(pheno_path, delimiter = "\t").loc[old_in_new, :]
PCs_info = pd.read_csv(PCs_path, delimiter = "\t").loc[old_in_new, :]
Betas_info = pd.read_csv(Betas_path, delimiter = "\t")

phenotypes_info_prev = pd.read_csv("QTL_phenotypes_eids.txt", delimiter = "\t", header = None)
if not np.all(phenotypes_info["eid"].to_numpy() == phenotypes_info_prev[0].to_numpy()):
    print("exiting: eids are not in the expected order")
    exit()
if not np.all(PCs_info["eid"].to_numpy() == phenotypes_info_prev[0].to_numpy()):
    print("exiting: eids are not in the expected order")
    exit()

phenotypes = phenotypes_info.loc[:, pheno].to_numpy()
Betas = Betas_info.loc[:, pheno].to_numpy()
PCs = PCs_info[PCs_info.columns[PCs_info.columns != "eid"]].to_numpy()
PCs2 = np.concatenate([PCs, np.ones((len(PCs), 1))], axis = 1)
LR_offset_all = np.matmul(PCs2, Betas)

path = "../step7_adjust_HF_for_covariates_logistic_PCA/X.txt"
gender_data = pd.read_csv(path, delimiter = "\t", usecols = ["eid", "22001-0.0"], dtype = int)
new_fam_df = pd.DataFrame(np.array([new_fam, np.arange(len(new_fam))]).T)
new_fam_df.columns = ["eid", "index"] 
gender_data = gender_data.merge(new_fam_df, on = "eid", how = "inner").sort_values(by = "index")
is_male = gender_data["22001-0.0"].to_numpy(dtype = "bool")


if chr == "Y":
    phenotypes = phenotypes[is_male]
    LR_offset_all = LR_offset_all[is_male]

prefix = "../step9_regress_phenotypes_against_SNPs_PCA/genotype_metadata/genotype_metadata_chr"
cols = ["SNP", "C(HOM A1)", "C(HET)", "C(HOM A2)", "C(MISSING)"]
all_missingness = pd.read_csv(prefix + chr + ".frqx", delimiter = "\t", usecols = ["SNP", "C(MISSING)"])
all_afs = pd.read_csv(prefix + chr + ".frq", delim_whitespace = True, usecols = ["SNP", "MAF"])
if not np.all(all_missingness["SNP"] == rsIDs) and np.all(all_afs["SNP"] == rsIDs):
    print("exiting: frq and frqx rsID are not ordered as expected")
all_missingness = all_missingness["C(MISSING)"].to_numpy()/len(new_fam_main)
all_afs = all_afs["MAF"].to_numpy()
all_mafs = maf2 = np.min(np.array([all_afs, 1 - all_afs]), axis = 0)
non_missing_indices = np.arange(len(rsIDs))[all_missingness <= 0.02]
non_constant_indices =  np.arange(len(rsIDs))[all_mafs >= 0.005]
good_indices = np.intersect1d(non_missing_indices, non_constant_indices)
intervals = np.cumsum([0] + int(len(good_indices)/50)*[50] + [len(good_indices)%50])

N = len(good_indices)
new_rsIDs = []
out = np.zeros((N, 3))
i_final = 0
for j in tqdm(range(j_start, len(intervals) - 1)):

    if j == j_start + 500:
        break

    col_indices = good_indices[intervals[j]:intervals[j + 1]]
    rsIDs_chunk = rsIDs[col_indices].to_numpy()
    new_rsIDs += rsIDs_chunk.tolist()

    # my_rsID = 'rs5747212'
    # if my_rsID in rsIDs_chunk:
    #     pdb.set_trace()
    #     genotypes = open_bed(geno_path  + ".bed", num_threads = 1).read(np.s_[sorted_indices, col_indices])
    #     ind = np.where(rsIDs_chunk == my_rsID)[0][0]
    #     zz = nanlinregress(genotypes[:, ind], phenotypes)
    # else:
    #     continue

    genotypes = open_bed(geno_path  + ".bed", num_threads = 1).read(np.s_[sorted_indices, col_indices])
    genotypes = genotypes[new_in_old, :]
    is_male2 = is_male
    if chr == "Y":
        genotypes, is_male2  = genotypes[is_male], is_male[is_male]

    afs = all_afs[col_indices]
    mafs = all_mafs[col_indices]
    af_is_not_maf = (afs != mafs)
    genotypes[:, af_is_not_maf] = 2 - genotypes[:, af_is_not_maf]
    valued_indices = np.logical_or(np.isnan(genotypes), np.isnan(phenotypes).reshape(-1,1)) == False

    if chr == "X":
        all_genotypes = [genotypes[is_male2], genotypes[is_male2 == False]]
        all_phenotypes = [phenotypes[is_male2], phenotypes[is_male2 == False]]
        LR_offset_all2 = [LR_offset_all[is_male2], LR_offset_all[is_male2 == False]]
        all_valued_indices = [valued_indices[is_male2], valued_indices[is_male2 == False]]
        all_vals = zip(all_genotypes, all_phenotypes, all_valued_indices, [0, 1], LR_offset_all2)
        L = len(genotypes[0])
        p_main, p_unadj = np.zeros((2, L)), np.zeros((2, L))
        num_groups = 2
    else:
        all_vals = zip([genotypes], [phenotypes], [valued_indices], [0], [LR_offset_all])
        num_groups = 1
    for geno, pheno_vals, val_inds, ind, LR_offset in all_vals:
        i = COPY(i_final)
        for k in range(len(geno[0])):

            g = geno[val_inds[:, k], k]
            offset = LR_offset[val_inds[:, k]]
            out[i, 2] = len(g)
            p = pheno_vals[val_inds[:, k]]
            g2_vals = (g == 2)
            is_standard = np.sum(g2_vals, axis = 0) > 1000 and chr != "Y"
            
            if not is_standard: 
                g[g2_vals] = 1

            g = np.array([g, np.ones(len(g))]).T
            family = sm.genmod.families.Binomial()
            family.link = sm.genmod.families.links.logit()
            model = sm.GLM(p, g, family = family, offset = offset)
            # This and model = sm.Logit(p, g) are equivalent to the 6th decimal
            # This is about 95% similar to LR with PCs as covariates: ~0.01 -> ~0.05
            # model = sm.Logit(p, g, offset = offset) does not include the PC offset
            p_val = model.fit(disp=0).pvalues[0]

            if num_groups == 1:
                out[i, 1] = p_val
            else:
                p_main[ind][k] = p_val
            i += 1

    if num_groups == 2:
        out[(i - k - 1):(i), 1] = chi2(4).sf(-2*np.sum(np.log(p_main), axis = 0))
    i_final = COPY(i)

o_df = pd.DataFrame(out[out[:, 2] != 0])
o_df.columns = ["rsID", "p_main", "N"]
o_df["rsID"] = new_rsIDs
suffix = "_chr" + chr + "_P" + str(pheno) + "_jstart" + str(j_start)
o_df.to_csv("binary_HF_QTL_output/binary_HF_QTL_effects" + suffix + ".txt", sep = "\t", header = True, index = False)