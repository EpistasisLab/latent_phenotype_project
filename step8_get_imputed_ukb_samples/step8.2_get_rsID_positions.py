import pandas as pd
import numpy as np
from tqdm import tqdm
import pdb

real_filtered_rsIDs_path = "../step4_remove_relatives/UKB_samples_unrelated.bim"
real_filtered_rsIDs = pd.read_csv(real_filtered_rsIDs_path, delim_whitespace = True, header = None)
real_filtered_rsIDs = real_filtered_rsIDs.loc[real_filtered_rsIDs[0] > 22, [1]]
real_filtered_rsIDs.to_csv("eid_filters/real_filtered_rsIDs.tab", sep = "\t", header = False, index = False)

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "XY"]
imputed_SNPs_path = "/project/UKB_moore/UKB_50978/imputed/penn_freeze_11132019/"
imputed_SNPs_file_names = [imputed_SNPs_path + "ukb_mfi_chr" + chr + "_v3.txt" for chr in chromosomes]

# Filtering low quality SNPs, SNPs with mafs less than 0.005, and SNPs that aren't bi-allelic
all_SNPs_chr = []
num_good_GWAS_SNPs = 0
for i, chr in tqdm(enumerate(chromosomes)):

    # column 7 is the info score, column 5 is the maf, and column 2 is the chromosome position
    try:
        SNP_info_scores = pd.read_csv(imputed_SNPs_file_names[i], delimiter = "\t", header = None)
    except:
        pdb.set_trace()
    SNP_info_scores.columns = ["name", 1, 2, 3, 4, 5, 6, 7]
    SNP_info_scores["name"] = SNP_info_scores.loc[:, "name"].str.split("_", expand = True)[0]

    # seperates non-biallelic SNPs
    duplicate_SNP_names = SNP_info_scores[SNP_info_scores.duplicated(["name"])]["name"]
    non_biallelic_SNPs1 = np.unique(SNP_info_scores[SNP_info_scores["name"].isin(duplicate_SNP_names)][1])
    duplicate_SNP_pos = SNP_info_scores[SNP_info_scores.duplicated([2])][2]
    non_biallelic_SNPs2 = np.unique(SNP_info_scores[SNP_info_scores[2].isin(duplicate_SNP_pos)][1])
    non_biallelic_SNPs3 = SNP_info_scores[SNP_info_scores.duplicated([1])][1].to_numpy()
    non_biallelic_SNPs = np.union1d(np.union1d(non_biallelic_SNPs1, non_biallelic_SNPs2), non_biallelic_SNPs3)
    all_SNPs_subset = SNP_info_scores[~SNP_info_scores[1].isin(non_biallelic_SNPs)]

    # I want at least an effective sample size (info score times sample size) of 50000 and a maf >= 0.005
    good_SNP_indices = np.logical_and(all_SNPs_subset[7] >= (50000/380940), all_SNPs_subset[5] >= 0.005)
    all_SNPs_chr.append(all_SNPs_subset[good_SNP_indices])

for i, SNP_set in enumerate(all_SNPs_chr): 
    out_file = "eid_filters/SNP_positions_chr" + chromosomes[i] + ".txt"
    SNP_set[1].to_csv(out_file, sep = "\t", header = False, index = False)