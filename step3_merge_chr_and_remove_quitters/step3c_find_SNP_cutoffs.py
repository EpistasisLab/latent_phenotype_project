import numpy as np
import pandas as pd
import os
import pdb
from matplotlib import pyplot as plt

hwe_path = "UKB_samples.hwe"
mis_path = "UKB_samples.lmiss"
maf_path = "UKB_samples.frq"
hwe_data = pd.read_csv(hwe_path, delim_whitespace = True, usecols = ["SNP", "P"], header = 0)
missingness_data = pd.read_csv(mis_path, delim_whitespace = True, usecols = ["SNP", "F_MISS"], header = 0)
maf_data = pd.read_csv(maf_path, delim_whitespace = True, usecols = ["SNP", "MAF"], header = 0)
all_data = maf_data.merge(hwe_data, on = "SNP")
data = all_data.merge(missingness_data, on = "SNP")
data = data[data["P"] > 0]
data["P"] = np.log10(data["P"])
hwe_logP_values = data["P"].to_numpy()
missingness = data["F_MISS"].to_numpy()
mafs = data["MAF"].to_numpy()

average_missingness = []
average_maf = []
average_missingness.append(np.mean(data[np.logical_and(data["MAF"] > 0.00001, data["MAF"] <= 0.001)]["F_MISS"]))
average_missingness.append(np.mean(data[np.logical_and(data["MAF"] > 0.001, data["MAF"] <= 0.01)]["F_MISS"]))
average_missingness.append(np.mean(data[np.logical_and(data["MAF"] > 0.01, data["MAF"] <= 0.1)]["F_MISS"]))
average_missingness.append(np.mean(data[np.logical_and(data["MAF"] > 0.1, data["MAF"] <= 0.5)]["F_MISS"]))
average_maf.append(np.mean(data[np.logical_and(data["MAF"] > 0.00001, data["MAF"] <= 0.001)]["MAF"]))
average_maf.append(np.mean(data[np.logical_and(data["MAF"] > 0.001, data["MAF"] <= 0.01)]["MAF"]))
average_maf.append(np.mean(data[np.logical_and(data["MAF"] > 0.01, data["MAF"] <= 0.1)]["MAF"]))
average_maf.append(np.mean(data[np.logical_and(data["MAF"] > 0.1, data["MAF"] <= 0.5)]["MAF"]))

low_missingness_SNPs = data[data["F_MISS"] < 0.02]
clean_data = low_missingness_SNPs[low_missingness_SNPs["P"] > -6]
clean_average_missingness = []
clean_average_maf = []
clean_average_missingness.append(np.mean(clean_data[np.logical_and(clean_data["MAF"] > 0.00001, clean_data["MAF"] <= 0.001)]["F_MISS"]))
clean_average_missingness.append(np.mean(clean_data[np.logical_and(clean_data["MAF"] > 0.001, clean_data["MAF"] <= 0.01)]["F_MISS"]))
clean_average_missingness.append(np.mean(clean_data[np.logical_and(clean_data["MAF"] > 0.01, clean_data["MAF"] <= 0.1)]["F_MISS"]))
clean_average_missingness.append(np.mean(clean_data[np.logical_and(clean_data["MAF"] > 0.1, clean_data["MAF"] <= 0.5)]["F_MISS"]))
clean_average_maf.append(np.mean(clean_data[np.logical_and(clean_data["MAF"] > 0.00001, clean_data["MAF"] <= 0.001)]["MAF"]))
clean_average_maf.append(np.mean(clean_data[np.logical_and(clean_data["MAF"] > 0.001, clean_data["MAF"] <= 0.01)]["MAF"]))
clean_average_maf.append(np.mean(clean_data[np.logical_and(clean_data["MAF"] > 0.01, clean_data["MAF"] <= 0.1)]["MAF"]))
clean_average_maf.append(np.mean(clean_data[np.logical_and(clean_data["MAF"] > 0.1, clean_data["MAF"] <= 0.5)]["MAF"]))

plt.plot(average_maf, average_missingness, 'x', label = "raw data") 
plt.plot(clean_average_maf, clean_average_missingness, 'o', label = "HWE log(p) > -6 and SNP missingness < 0.02")
plt.xlabel("average maf", fontsize = 12) 
plt.ylabel("average missingness", fontsize = 12)
plt.legend()
plt.savefig("maf_vs_missingness.png")
plt.clf()

plt.hist(missingness, bins = 500)
plt.xlabel("missingness", fontsize = 12)
plt.xlim(0, 0.2)
plt.ylabel("frequency", fontsize = 12)
plt.savefig("missingness_distribution.png")
plt.clf()

plt.hist(hwe_logP_values, bins = 500)
plt.xlabel("hwe_logP_values", fontsize = 12)
plt.xlim(-20, 0)
plt.ylabel("frequency", fontsize = 12)
plt.savefig("hwe_logP_values_distribution.png")
plt.clf()
