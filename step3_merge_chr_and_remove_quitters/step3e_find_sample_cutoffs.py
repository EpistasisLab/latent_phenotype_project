import numpy as np
import pandas as pd
import pdb
import os
from copy import deepcopy as COPY
from matplotlib import pyplot as plt

# imports the data
het_path = "UKB_samples_half_filtered.het"
mis_path = "UKB_samples_half_filtered.imiss"
sex_path = "UKB_samples_half_filtered.sexcheck"
heterozygosity_data = pd.read_csv(het_path, delim_whitespace=True, usecols = ["FID", "O(HOM)", "N(NM)"], header = 0)
missingness_data = pd.read_csv(mis_path, delim_whitespace=True, usecols = ["FID", "F_MISS"], header = 0)
sex_data = pd.read_csv(sex_path, delim_whitespace=True, usecols = ["FID", "PEDSEX", "SNPSEX", "STATUS", "F"], header = 0)
data = heterozygosity_data.merge(missingness_data, on = "FID")
data = data.merge(sex_data, on = "FID")
data["heterozygosity"] = ((data["N(NM)"] - data["O(HOM)"])/data["N(NM)"]).to_numpy()

# determines heterozygosity outlier bounds
het_vals = data["heterozygosity"].to_numpy()
low_het_cutoff = np.mean(het_vals) - 3*np.std(het_vals)
high_het_cutoff = np.mean(het_vals) + 3*np.std(het_vals)
plt.hist(het_vals, bins = 1000, label = "heterozygosity distribution cutoffs")
plt.axvline(low_het_cutoff)
plt.axvline(high_het_cutoff)
plt.legend()
plt.xlabel("heterozygosity rate", fontsize = 12)
plt.ylabel("number of value instances", fontsize = 12)
plt.savefig("heterozygosity.png")
plt.clf()

# determines gender heterozygosity outlier bounds
males, females = data[data["PEDSEX"] == 1], data[data["PEDSEX"] == 2]
male_F_scores, female_F_scores = males["F"], females["F"]
# females have a distribution of X heterozygosity, so outliers need to be determined. 
# males should have no X heterozygosity (i.e. they have a Y chr), so outliership is just F < 1 
low_F_cutoff = np.mean(female_F_scores) - 3*np.std(female_F_scores)
high_F_cutoff = np.mean(female_F_scores) + 3*np.std(female_F_scores)
plt.hist(female_F_scores, bins = 1000, label = "female F score distribution cutoffs")
plt.axvline(low_F_cutoff)
plt.axvline(high_F_cutoff)
plt.legend()
plt.xlabel("gender heterozygosity F score", fontsize = 12)
plt.ylabel("number of value instances", fontsize = 12)
plt.savefig("heterozygosity_X_chr.png")
plt.clf()
 
# removes all unwanted samples
het_outlier_indices = np.logical_or(het_vals < low_het_cutoff, het_vals > high_het_cutoff)
het_outlier_FIDs  = data.loc[het_outlier_indices, "FID"].to_numpy()
miss_outlier_indices = data["F_MISS"] > 0.02
miss_outlier_FIDs = data.loc[miss_outlier_indices, "FID"].to_numpy()
female_outlier_indices = np.logical_or(female_F_scores < low_F_cutoff, female_F_scores > high_F_cutoff)
female_outlier_FIDs = females.loc[female_outlier_indices , "FID"].to_numpy()
male_outlier_indices = males["F"] != 1
male_outlier_FIDs = males.loc[male_outlier_indices , "FID"].to_numpy()
low_quality_FIDs1 = np.union1d(het_outlier_FIDs, miss_outlier_FIDs)
low_quality_FIDs2 = np.union1d(female_outlier_FIDs, male_outlier_FIDs)
low_quality_FIDs = np.union1d(low_quality_FIDs1, low_quality_FIDs2)

#------------------------------------------------------------------------------------------------------------------------------------
# start of step 3 methodology change for reviewers
#------------------------------------------------------------------------------------------------------------------------------------
non_WB_eids = pd.read_csv("non_white_british_eids.tab", sep = "\t", header = None)[0].to_numpy()
low_quality_FIDs = np.union1d(low_quality_FIDs, non_WB_eids)
#------------------------------------------------------------------------------------------------------------------------------------
# end of step 3 methodology change for reviewers
#------------------------------------------------------------------------------------------------------------------------------------

low_quality_FIDs = pd.DataFrame(np.array([low_quality_FIDs, low_quality_FIDs]).T)
low_quality_FIDs.to_csv("low_quality_FIDs.txt", sep = "\t", header = False, index = False)

