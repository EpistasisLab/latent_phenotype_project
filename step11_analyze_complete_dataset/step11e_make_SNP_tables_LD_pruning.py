import numpy as np
import pandas as pd
import os
from bed_reader import open_bed
from functools import reduce
import pdb

output_names = ["snps_output/" + name for name in os.listdir("snps_output")]
output_files = [pd.read_csv(name, delimiter = "\t") for name in output_names]
output_types = np.array([file.columns[0] for file in output_files])

# I've manually selected a subset of terms with potentially AHF associated word fragments
# they are based on manual examination of the "GWAS traits" column in "output_files" dataframes
possible_AHF_terms = ["heart", "cardi", "vessel", "vein", "arter"]
possible_AHF_terms += ["capillar", "aort", "aneurysm", "calc", "clot"]
possible_AHF_terms += ["cholest", "lipopro", "ldl", "hdl", "vldl"]
possible_AHF_terms += ["diabet", "hypertens", "syst", "diast", "pressure"]
possible_AHF_terms += ["arrhythmia", "atheroscl", "insulin", "stroke", "embol"]
possible_AHF_terms += ["thromb", "coagul", "leptin", "ischem", "platelet"]
possible_AHF_terms += ["renin", "angiotensin", "sphingomyelin", "triacylglyceride", "venous"]

no_info_inds = np.array([len(info) for info in output_files]) == 0
SNPs_with_info = pd.concat(np.array(output_files, dtype = object)[no_info_inds == False])
unique_traits = np.unique([i.lower() for i in SNPs_with_info["GWAS Trait"].tolist()]).astype(str)
term_sets = [[trait for trait in unique_traits if term in trait] for term in possible_AHF_terms]
term_sets = reduce(np.union1d, term_sets)
output_files = [df if len(df) == 0 else df.loc[df["GWAS Trait"].str.lower().isin(term_sets), :] for df in output_files]

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
 
# old formatting
# counts as opposed to SNP_counts count the number of SNP effects as opposed to the number of unique SNPs.
# counts = pd.DataFrame([categories, novel_counts, novel_SNP_counts, known_counts, known_SNP_counts, invalid_counts]).T

counts = pd.DataFrame([categories, novel_SNP_counts, known_SNP_counts]).T
counts = counts.sort_values(0, ascending = False)
counts.columns = ["effect type", "novel SNPs", "known SNPs"]

counts_logistic_PCA = counts[["novel SNPs", "known SNPs"]].to_numpy()[0:4].reshape(-1)
counts_PCA = counts[["novel SNPs", "known SNPs"]].to_numpy()[4:8].reshape(-1)
counts_NN = counts[["novel SNPs", "known SNPs"]].to_numpy()[8:12].reshape(-1)
counts2 = pd.DataFrame([counts_PCA, counts_logistic_PCA, counts_NN])
counts2[8] = ["PCA", "logistic PCA", "autoencoder"]
counts2 = counts2[[8, 0, 1, 2, 3, 4, 5, 6, 7]]
cols = ["latent phenotype model", "Main_Novel_SNP", "Main_Known_SNP"]
cols += ["GxSmoking_Novel_SNP", "GxSmoking_Known_SNP"]
cols += ["GxGender_Novel_SNP", "GxGender_Known_SNP"]
cols += ["GxAlcohol_Novel_SNP", "GxAlcohol_Known_SNP"]
counts2.columns = cols
counts2.to_csv("table1b.txt", sep = "\t", header = True, index = False)
