import numpy as np
import pandas as pd
import os
from functools import reduce
from scipy.stats import chi2
from copy import deepcopy as COPY
import pdb


base = "../step9_regress_phenotypes_against_SNPs_PCA/"
num_pheno = 16
pg = (5E-8)/16
chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
env_factors = ["_gender", "_exercise", "_smoking", "_alcohol"]

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

        try:
            significant_QTLs = pd.concat(significant_QTLs)
        except: 
            pdb.set_trace()
        sig_QTL_path = "significant_QTLs" + name + "/significant_QTLs_P" + str(i) + ".txt"
        all_significant_QTLs.append(significant_QTLs)
        significant_QTLs.to_csv(sig_QTL_path, sep = "\t", header = True, index = False)

    all_significant_rsIDs = pd.DataFrame(np.unique(pd.concat(all_significant_QTLs)["rsID"]))
    all_significant_rsIDs.to_csv("significant_rsIDs" + name + ".txt", sep = "\t", header = False, index = False)
