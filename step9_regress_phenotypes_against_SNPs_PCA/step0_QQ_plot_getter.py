import numpy as np
import pandas as pd
import argparse
import os
import pdb
from scipy.stats import beta
from scipy.stats import chi2
from tqdm import tqdm
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--env', nargs = 1, type = str, action = "store", dest = "env")
parser.add_argument('--pheno', nargs = 1, type = str, action = "store", dest = "pheno")
parser.add_argument('--GWAS', nargs = '*', type = str, action = "store", dest = "GWAS")
args = parser.parse_args()
env_name = (args.env)[0]
pheno = (args.pheno)[0]
GWAS = args.GWAS

if (not os.path.isdir("QQ_plots_" + env_name)) and pheno == '0':
    os.mkdir("QQ_plots_" + env_name)

chrs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"]
chrs += ["15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]

if GWAS == ["normal"]:
    prefix = "binary_HF_QTL_output"
    file_names = [prefix + "/" + n for n in os.listdir(prefix) if pheno in n]
    file = pd.concat([pd.read_csv(n, delimiter = "\t") for n in file_names])
else:
    prefix = "QTL_output_" + env_name + "/QTL_effects_chr"
    file_names = [prefix + j + "_P" + pheno + ".txt" for j in chrs]
    file = pd.concat([pd.read_csv(n, delimiter = "\t") for n in file_names])

N = len(file)
medians = beta.ppf(0.5, np.arange(1, N+1), N - np.arange(N))
medians = -np.log10(np.flip(medians))
ubs = beta.ppf(0.975, np.arange(1, N+1), N - np.arange(N))
ubs = -np.log10(np.flip(ubs))
lbs = beta.ppf(0.025, np.arange(1, N+1), N - np.arange(N))
lbs = -np.log10(np.flip(lbs))

if GWAS == ["normal"]:
    pvals = -np.log10(np.flip(np.sort(file["p_main"].to_numpy())))
else:
    pvals = -np.log10(np.flip(np.sort(file["p_null2"].to_numpy())))


chi_alt = chi2.isf(10**(-np.median(pvals)), 1, loc=0, scale=1)
chi_null = chi2.isf(10**(-np.median(medians)), 1, loc=0, scale=1)
inflation_factor = pd.DataFrame([chi_alt/chi_null])
name = "QQ_plots_" + env_name + "/phenotype" + str(pheno) + "_IF.txt"
inflation_factor.to_csv(name, sep = "\t", header = False, index = False)

fig, ax = plt.subplots(figsize = (10, 7.5))
plt.title('p-value inflation curves', fontsize=20)

sorted_p_vals = np.sort(pvals)
pvals_cdf = np.cumsum(np.ones(len(sorted_p_vals))/len(sorted_p_vals))
ax.plot(sorted_p_vals, pvals_cdf, "b-")
ax.plot(sorted_p_vals, np.ones(len(pvals_cdf)), "k-")

ax2 = ax.twinx()
ax2.plot(medians, pvals, "g-")
ax2.plot(medians, lbs, "r-")
ax2.plot(medians, ubs, "r-")
ax2.plot(medians, medians, "k-")
ax2.fill_between(medians, lbs, ubs, color = "r", alpha = 0.5)

ax.set_xlim([0, np.max(medians)])
ax.set_xlabel("expected p values under null hypothesis (95% CI)", color = 'r', fontsize=20)
ax.set_ylabel("CDF of actual p-values", color = 'b', fontsize=20)
ax2.set_ylabel("actual p-values (QQ plot)", color = 'g', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize = 20)
ax2.tick_params(axis='both', which='major', labelsize = 20)
plt.tight_layout()
plt.savefig("QQ_plots_" + env_name + "/QQ_plot_pheno" + pheno + ".png")
plt.clf()