import numpy as np
import pandas as pd
from bed_reader import open_bed
import os
import pdb

output_names = ["snps_output/" + name for name in os.listdir("snps_output")]
output_files = [pd.read_csv(name, delimiter = "\t") for name in output_names]
output_types = np.array([file.columns[0] for file in output_files])

input_rsIDs = np.array([name.split("_")[3] for name in output_names])
input_E = np.array([name.split("/")[1].split("_")[0] for name in output_names])
lengths = np.array([len(file) for file in output_files])

valid_output = np.array(output_files, dtype = object)
invalid_output = valid_output[output_types == 'invalid name']
valid_output = valid_output[output_types != 'invalid name']
invalid_rsIDs = input_rsIDs[output_types == 'invalid name']
valid_rsIDs = input_rsIDs[output_types != 'invalid name']
invalid_E = input_E[output_types == 'invalid name']
valid_E = input_E[output_types != 'invalid name']
invalid_lengths = lengths[output_types == 'invalid name']
valid_lengths = lengths[output_types != 'invalid name']
output_types = output_types[output_types != 'invalid name']

paths = ["all_sig_rsIDs_PCA", "all_sig_rsIDs_logistic_PCA", "all_sig_rsIDs_NN"]


novel_counts = []
novel_SNP_counts = []
known_counts = []
known_SNP_counts = []
invalid_counts = []
categories = []
for E in np.unique(input_E):

    invalid_rsIDs_E = invalid_rsIDs[invalid_E == E]
    invalid_rsIDs_E = np.array([i.replace("!", ":") for i in invalid_rsIDs_E])
    
    valid_rsIDs_E = valid_rsIDs[valid_E == E]
    valid_lengths_E = valid_lengths[valid_E == E]
    known_rsIDs = valid_rsIDs_E[valid_lengths_E > 0]
    novel_rsIDs = valid_rsIDs_E[valid_lengths_E == 0]

    valid_output_E = valid_output[valid_E == E]
    known_rsID_traits = valid_output_E[valid_lengths_E > 0]

    for path in paths:
        
        fname = path + "/rsIDs_" + E + "_effects.txt"
        snp_set = pd.read_csv(fname, delimiter = "\t")
        snp_set_rsIDs = snp_set["rsID"].to_numpy()
        snp_set_rsIDs = np.array([i.split("_")[0] for i in snp_set_rsIDs])
        snp_set["rsID"] = snp_set_rsIDs 
        
        snp_subset_known = pd.DataFrame(known_rsIDs, columns = ["rsID"])
        snp_subset_known = snp_subset_known.merge(snp_set, on = "rsID", how = "inner")
        snp_subset_novel = pd.DataFrame(novel_rsIDs, columns = ["rsID"])
        snp_subset_novel = snp_subset_novel.merge(snp_set, on = "rsID", how = "inner")
        snp_subset_invalid = pd.DataFrame(invalid_rsIDs_E, columns = ["rsID"])
        snp_subset_invalid = snp_subset_invalid.merge(snp_set, on = "rsID", how = "inner")

        known_counts.append(len(snp_subset_known))
        novel_counts.append(len(snp_subset_novel))
        invalid_counts.append(len(snp_subset_invalid))
        categories.append(path.split("rsIDs_")[1] + "_" + E)

        out_fname1 = "step11e_" + path.split("rsIDs_")[1] + "_" + E + "_rsIDs_known.txt"
        out_fname2 = "step11e_" + path.split("rsIDs_")[1] + "_" + E + "_rsIDs_novel.txt"
        df1 = snp_subset_known["rsID"].drop_duplicates()
        df2 = snp_subset_novel["rsID"].drop_duplicates()

        plink_prefix = "../step10_get_significant_SNPs_" + path.split("rsIDs_")[1] + "/"
        if E == "GxSmoking" or E == "main":
            plink_path = plink_prefix + "significant_SNPs_plink_files_smoking"
        if E == "GxAlcohol":
            plink_path = plink_prefix + "significant_SNPs_plink_files_alcohol"
        if E == "GxGender" or E == "main":
            plink_path = plink_prefix + "significant_SNPs_plink_files_gender"
        bed_paths = [plink_path + "/" + i for i in os.listdir(plink_path) if ".bed" in i]
        bim_paths = [plink_path + "/" + i for i in os.listdir(plink_path) if ".bim" in i]
        bim_files = pd.concat([pd.read_csv(path, delimiter = "\t", header = None) for path in bim_paths])

        known_SNP_count = 0
        known_independent_rsIDs = []
        novel_SNP_count = 0
        novel_independent_rsIDs = []
        for bim_path, bed_path in zip(bim_paths, bed_paths):
            bim_file = pd.read_csv(bim_path, delimiter = "\t", header = None)
            bim_inds_known = bim_file[1].isin(df1.to_numpy()).to_numpy()
            bim_inds_novel = bim_file[1].isin(df2.to_numpy()).to_numpy()
            for l, inds in enumerate([bim_inds_known, bim_inds_novel]): 
                ind_sum = np.sum(inds)
                rsIDs = bim_file[1].to_numpy()[inds]
                if ind_sum < 2:
                    if l == 0: 
                        known_SNP_count += ind_sum
                        known_independent_rsIDs += rsIDs.tolist()
                    if l == 1: 
                        novel_SNP_count += ind_sum
                        novel_independent_rsIDs += rsIDs.tolist()
                else:
                    X = open_bed(bed_path).read(np.s_[:, inds])
                    for i in range(ind_sum): X[np.isnan(X[:, i]), i] = np.nanmean(X[:, i])
                    updates = True
                    while updates == True:
                        max_corrs = np.max((np.corrcoef(X.T) - np.eye(len(X[0])))**2, axis = 0)
                        max_corr = np.max(max_corrs)
                        if max_corr < 0.8:
                            updates = False
                            if l == 0: 
                                known_SNP_count += len(X[0])
                                known_independent_rsIDs += rsIDs.tolist()
                            if l == 1: 
                                novel_SNP_count += len(X[0])
                                novel_independent_rsIDs += rsIDs.tolist()
                        else:
                            rsIDs = rsIDs[np.setdiff1d(range(len(X[0])), np.argmax(max_corrs))]
                            X = X[:, np.setdiff1d(range(len(X[0])), np.argmax(max_corrs))]

        novel_SNP_counts.append(novel_SNP_count)
        known_SNP_counts.append(known_SNP_count)
        known_df = pd.DataFrame(known_independent_rsIDs, columns = ["rsID"])
        novel_df = pd.DataFrame(novel_independent_rsIDs, columns = ["rsID"])
        known_df.to_csv(out_fname1, sep = "\t", header = True, index = False)
        novel_df.to_csv(out_fname2, sep = "\t", header = True, index = False)
 

counts = pd.DataFrame([categories, novel_counts, novel_SNP_counts, known_counts, known_SNP_counts, invalid_counts]).T
counts = counts.sort_values(0)
counts.columns = ["effect type", "novel SNP hits", "novel SNPs", "known SNP hits", "known SNPs", "invalid names"]
counts.to_csv("step11e_counts.txt", sep = "\t", header = True, index = False)
