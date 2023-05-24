import numpy as np
import pandas as pd
import os
from functools import reduce
from scipy.stats import chi2
from copy import deepcopy as COPY
from matplotlib import pyplot as plt
import pdb

base = "../step9_regress_phenotypes_against_SNPs_logistic_PCA/"
pg = (5E-8)/16
num_pheno = 16
chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
env_factors = ["_alcohol", "_exercise", "_gender", "_smoking"]

phenotypes = pd.read_csv(base + "QTL_phenotypes.txt", delimiter = "\t", header = None) 
for i in range(num_pheno):
    plt.hist(phenotypes.to_numpy()[:, i], bins = 1000)
    plt.savefig(base + "phenotype" + str(i) + ".png")
    plt.clf()
    plt.hist(phenotypes.to_numpy()[:, i], bins = 1000)
    plt.ylim((-10, 10000))
    plt.savefig(base + "phenotype" + str(i) + "_truncated.png")
    plt.clf()

for name in env_factors:

    if not os.path.exists("significant_QTLs" + name):
        os.mkdir("significant_QTLs" + name)

    all_significant_QTLs = []
    for i in range(num_pheno):
        significant_QTLs = []
        for c in chromosomes:

            QTL_path = base + "QTL_output" + name + "/QTL_effects_chr" + c + "_P" + str(i) + ".txt"
            try:
                QTLs = pd.read_csv(QTL_path, delimiter = "\t", header = 0)
            except:
                continue
            QTLs["chr"] = c
            QTL_is_significant = QTLs["p_main"].to_numpy() < pg
            significant_QTLs.append(QTLs.loc[QTL_is_significant, :])

        significant_QTLs = pd.concat(significant_QTLs)
        sig_QTL_path = "significant_QTLs" + name + "/significant_QTLs_P" + str(i) + ".txt"
        all_significant_QTLs.append(significant_QTLs)
        significant_QTLs.to_csv(sig_QTL_path, sep = "\t", header = True, index = False)

    all_significant_rsIDs = pd.DataFrame(np.unique(pd.concat(all_significant_QTLs)["rsID"]))
    all_significant_rsIDs.to_csv("significant_rsIDs" + name + ".txt", sep = "\t", header = False, index = False)
