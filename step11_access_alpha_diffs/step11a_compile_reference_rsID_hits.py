import numpy as np 
import pandas as pd
import tabula
import pdb
import requests
import sys
import os
from copy import deepcopy as COPY
from bed_reader import open_bed
from scipy.stats import chi2
from scipy.stats import pearsonr
from scipy.stats import fisher_exact
from scipy.stats import mannwhitneyu

# prereqs:
# conda install -c cyclus java-jre
# pip install tabula-py

chromosomes = list(np.arange(1, 23).astype(str)) + ["X", "XY", "Y", "MT"]
bim_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
bim_paths = [bim_prefix + i + ".bim" for i in chromosomes]
bim_files = [pd.read_csv(path, delim_whitespace = True, header = None) for path in bim_paths]
bim_files[-4] = pd.concat([bim_files[-4], bim_files[-3]]).sort_values(by = 3)
bim_files[-4].loc[:, 0] = "X"

my_hits_path_male = "../step10_get_significant_SNPs/QTL_hits_male.txt"
my_hits_path_female = "../step10_get_significant_SNPs/QTL_hits_female.txt"
my_hits_male = pd.read_csv(my_hits_path_male, delimiter = "\t", usecols = ['rsID', 'p', 'alpha', 'phenotype', "chr"])
my_hits_male = my_hits_male.loc[my_hits_male["phenotype"] == 0, :]
del my_hits_male["phenotype"]
my_hits_male.columns = ['rsID', 'p_m', 'alpha_m', "chr"]
my_hits_female = pd.read_csv(my_hits_path_female, delimiter = "\t", usecols = ['rsID', 'p', 'alpha', 'phenotype', "chr"])
my_hits_female = my_hits_female.loc[my_hits_female["phenotype"] == 0, :]
del my_hits_female["phenotype"]
my_hits_female.columns = ['rsID', 'p_f', 'alpha_f', "chr"]

p_vals_male = my_hits_male["p_m"].to_numpy()
p_vals_female = my_hits_female["p_f"].to_numpy()
chi2_vals = -2*(np.log(p_vals_male) + np.log(p_vals_female))
p_vals_joint = chi2.sf(chi2_vals, df = 4)
my_hits_all = my_hits_male.merge(my_hits_female, on = ["rsID", "chr"], how = "inner")
my_hits_all.loc[:, "p_joint"] = p_vals_joint

my_hits_chr = [my_hits_all.loc[my_hits_all["chr"] == i, :] for i in chromosomes]
my_hits_chr[-4] = pd.concat([my_hits_chr[-4], my_hits_chr[-3]])
my_hits_chr[-4].loc[:, "chr"] = "X"
void = my_hits_chr.pop(-3), bim_files.pop(-3), chromosomes.pop(-3)

chromosomes = list(np.arange(1, 23).astype(str)) + ["X", "Y", "MT"]
map_paths = ["step0/chr" + i + "_mapper.txt" for i in chromosomes]
mappers = [pd.read_csv(path, delimiter = "\t", header = None, dtype = int) for path in map_paths]

'''
# LATENT PHENOTYPE 0

# SNPs to consider (note: NOT does not mean against. Against is represented by a negative coefficient):
# Essential (primary) hypertension (I10, 0.9617654580108)
# Chronic ischaemic heart disease, unspecified (I259, 0.1034800682532) (NOT atherosclerotic)
# Other forms of chronic ischaemic heart disease (I258, 0.0991737521949) (NOT atherosclerotic) [no hits]
# Angina pectoris, unspecified (I209, 0.0894687878148) (NOT unstable angina, meaning NOT atherosclerotic) [no study]
# Hypertensive renal disease with renal failure (I120, 0.0775213949367) [no study]
# Atrial fibrillation and flutter (I48, 0.0688964032082)
# Peripheral vascular disease, unspecified (I739, 0.0523466869538)

# ---------------------- not included in SNP analysis
# Atherosclerotic heart disease (I251, 0.0487479916636)
# (Any_HF, 0.0454336694279)
# Stroke, not specified as haemorrhage or infarction (I64, 0.0326380558888)
# ---------------------- not included in SNP analysis

# LATENT ENV FACTOR 0

# Average night-time sound level of noise pollution  (24022, 0.7345352534766)
# Average daytime sound level of noise pollution  (24020, -0.671499153313)
# Average evening sound level of noise pollution  (24021, -0.064040225723)
# On blood pressure medication (6177 == 2, 0.0602221353211)
# On cholesterol lowering medication (6177 == 1, 0.0196378023182)

'''

df_names = ["I10.tsv", "I48.tsv", "I259.tsv", "I739.tsv", "24022.tsv", "24021.tsv", "24020.tsv"]
dfs = [pd.read_csv(name, delimiter = "\t") for name in df_names]
sig_SNPs = pd.concat(dfs).drop_duplicates(subset = "Lead Variant")["Lead Variant"]
chr_pos_a1_a2 = sig_SNPs.str.split("_", expand = True)
chr_pos_a1_a2_vec= [chr_pos_a1_a2[chr_pos_a1_a2[0] == str(i)] for i in chromosomes]

hit_rsIDs = []
hg37_misses = []
bim_misses = []
SNP_comparisons = []
geno_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output"
for mapper, cp, bim, chr, my_hits in zip(mappers, chr_pos_a1_a2_vec, bim_files, chromosomes, my_hits_chr):
    if len(cp) > 0 and len(my_hits) > 0:

        hg37_lbs, hg37_ubs, hg38_lbs, hg38_ubs, diffs = mapper.to_numpy().T
        chr, old_positions = cp[0].to_numpy()[0], cp[1].to_numpy(dtype = int)
        lb_indices = np.searchsorted(hg38_lbs, old_positions, side = "right") - 1
        ub_indices = np.searchsorted(hg38_ubs, old_positions, side = "left")
        if np.any(lb_indices != ub_indices):
            hg37_misses.append(cp.to_numpy()[lb_indices != ub_indices])
        new_positions = (old_positions - diffs[ub_indices])[lb_indices == ub_indices]

        rsID_bim_indices = np.where(bim[3].isin(new_positions))[0]
        hit_rsID_set = bim.loc[rsID_bim_indices, 1].to_numpy()
        if len(lb_indices) != len(rsID_bim_indices) + np.sum(lb_indices != ub_indices): 
            bim_misses.append(cp.to_numpy()[np.isin(new_positions, bim[3]) == False])
        hit_rsIDs.append(hit_rsID_set)

        my_hit_bim_indices = np.where(bim[1].isin(my_hits["rsID"]))[0]
        my_rsID_set = bim.loc[my_hit_bim_indices, 1].to_numpy()
        geno_path = geno_path_prefix + "/UKB_samples_chr" + chr
        geno_hits = open_bed(geno_path  + ".bed", num_threads = 1).read(np.s_[:, rsID_bim_indices])
        geno_hits_mine = open_bed(geno_path  + ".bed", num_threads = 1).read(np.s_[:, my_hit_bim_indices])
        max_corrs = []
        SNP_comparison = []
        for g_mine, rs_mine in zip(geno_hits_mine.T, my_rsID_set):
            all_corrs = []
            all_p = []
            for g in geno_hits.T:
                val_indices = np.logical_or(np.isnan(g), np.isnan(g_mine)) == False
                corr, p = pearsonr(g[val_indices], g_mine[val_indices])
                all_corrs.append(corr**2)
                all_p.append(p)
            max_r2_ind = np.argmax(np.array(all_corrs)**2)
            comparison = [chr, rs_mine, hit_rsID_set[max_r2_ind], 
                          all_corrs[max_r2_ind]**2, all_p[max_r2_ind]]
            SNP_comparison.append(comparison)
        SNP_comparison = pd.DataFrame(SNP_comparison)
        SNP_comparison.columns = ["chr", "rsID", "highest_LD_hit", "r_sq", "p"] 
        SNP_comparisons.append(SNP_comparison)

pdb.set_trace()
all_comparisons = pd.concat(SNP_comparisons)
novel_SNPs = (all_comparisons.loc[all_comparisons["r_sq"] < 0.1, :])
known_SNPs = (all_comparisons.loc[all_comparisons["r_sq"] >= 0.1, :])
novel_SNPs_info = novel_SNPs.merge(my_hits_all, on = "rsID", how = "inner")
known_SNPs_info = known_SNPs.merge(my_hits_all, on = "rsID", how = "inner")
novel_SNPs_info.sort_values(by = "r_sq").to_csv("novel_SNPs.txt", sep = "\t", index = False)
known_SNPs_info.sort_values(by = "r_sq").to_csv("known_SNPs.txt", sep = "\t", index = False)

# alphas exactly equal to -1 are removed because that is what I assign when no 2 genotypes are present. 
known_keepers = np.logical_and(known_SNPs_info["alpha_m"] != -1, known_SNPs_info["alpha_f"] != -1) 
novel_keepers = np.logical_and(novel_SNPs_info["alpha_m"] != -1, novel_SNPs_info["alpha_f"] != -1)
known_test_data = known_SNPs_info.loc[known_keepers, :]
novel_test_data = novel_SNPs_info.loc[novel_keepers, :]
log_p_ratio_known = (np.log10(known_test_data["p_m"])/np.log10(known_test_data["p_f"])).to_numpy()
log_p_ratio_known = np.max([log_p_ratio_known, 1/log_p_ratio_known], axis = 0)
log_p_ratio_novel = (np.log10(novel_test_data["p_m"])/np.log10(novel_test_data["p_f"])).to_numpy()
log_p_ratio_novel = np.max([log_p_ratio_novel, 1/log_p_ratio_novel], axis = 0)
known_test_data = known_test_data[log_p_ratio_known  < 3]
novel_test_data = novel_test_data[log_p_ratio_novel  < 3]
known_alpha_diffs = np.abs(known_test_data["alpha_m"] - known_test_data["alpha_f"])
novel_alpha_diffs = np.abs(novel_test_data["alpha_m"] - novel_test_data["alpha_f"])
mannwhitneyu(novel_alpha_diffs, known_alpha_diffs, alternative = "greater")
row1 = [N_known_alt_alpha, len(known_SNPs_info) - N_known_alt_alpha]
row2 = [N_novel_alt_alpha, len(novel_SNPs_info) - N_novel_alt_alpha]
fisher_exact([row1, row2])
print(1)