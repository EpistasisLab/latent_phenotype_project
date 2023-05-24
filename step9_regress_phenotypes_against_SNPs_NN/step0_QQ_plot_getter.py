import numpy as np
import pandas as pd
import argparse
import os
import pdb
from scipy.stats import beta
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

plt.plot(medians, pvals, "k-", label = "QQ plot")
plt.plot(medians, lbs, "r-", label = "expected 2.5th percentile")
plt.plot(medians, ubs, "b-", label = "expected 97.5th percentile")
plt.plot(medians, medians, "m-", label = "expected 50th percentile")
plt.fill_between(medians, lbs, ubs, alpha = 0.5)
plt.legend()
plt.xlabel("expected quantile p values")
plt.xlabel("actual quantile p values")
if not os.path.isdir("QQ_plots_" + env_name):
    os.mkdir("QQ_plots_" + env_name)
plt.savefig("QQ_plots_" + env_name + "/QQ_plot_pheno" + pheno + ".png")
plt.clf()