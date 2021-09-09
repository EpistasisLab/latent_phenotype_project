import numpy as np
import pandas as pd
import os
from functools import reduce
import pdb

if not os.path.exists("significant_QTLs"):
    os.mkdir("significant_QTLs")
if not os.path.exists("significant_vQTLs"):
    os.mkdir("significant_vQTLs")
base = "../step9_regress_phenotypes_against_SNPs/"
pg = (5E-8)/32

chromosomes = ["1", "2", "3", "4", "5", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]

all_significant_rsIDs = []
for i in range(16):
    significant_QTLs = []
    significant_vQTLs = []
    for c in chromosomes:
        QTL_path = base + "QTL_output/SNP_effects_chr" + c + "_P" + str(i) + ".txt"
        QTLs = pd.read_csv(QTL_path, delimiter = "\t", header = 0)
        QTLs["chr"] = c
        sig_QTLs = QTLs.loc[QTLs["p"] < pg, :]
        significant_QTLs.append(sig_QTLs)
        all_significant_rsIDs.append(sig_QTLs["rsID"].to_numpy())

        vQTL_path = base + "vQTL_output/vQTL_effects_chr" + c + "_P" + str(i) + ".txt"
        vQTLs = pd.read_csv(vQTL_path, delimiter = "\t", header = 0)
        vQTLs["chr"] = c
        sig_vQTLs = vQTLs.loc[vQTLs["p"] < pg, :]
        significant_vQTLs.append(sig_vQTLs)
        all_significant_rsIDs.append(sig_vQTLs["rsID"].to_numpy())

    significant_QTLs = pd.concat(significant_QTLs)
    sig_QTL_path = "significant_QTLs/significant_QTLs_P" + str(i) + ".txt"
    significant_QTLs.to_csv(sig_QTL_path, sep = "\t", header = True, index = False)

    significant_vQTLs = pd.concat(significant_vQTLs)
    sig_vQTL_path = "significant_vQTLs/significant_vQTLs_P" + str(i) + ".txt"
    significant_vQTLs.to_csv(sig_vQTL_path, sep = "\t", header = True, index = False)

all_significant_rsIDs = pd.DataFrame(reduce(np.union1d, all_significant_rsIDs))
all_significant_rsIDs.to_csv("significant_rsIDs.txt", sep = "\t", header = False, index = False)

