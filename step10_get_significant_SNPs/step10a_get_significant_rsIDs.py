import numpy as np
import pandas as pd
import os
from functools import reduce
from scipy.stats import chi2
from copy import deepcopy as COPY
import pdb

if not os.path.exists("significant_QTLs"):
    os.mkdir("significant_QTLs")
if not os.path.exists("significant_vQTLs"):
    os.mkdir("significant_vQTLs")
if not os.path.exists("GxE_candidate_rsIDs"):
    os.mkdir("GxE_candidate_rsIDs")
base = "../step9_regress_phenotypes_against_SNPs/"
pg = (5E-8)/32

chromosomes = ["1", "2", "3", "4", "5", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]

all_significant_rsIDs_male = []
all_significant_rsIDs_female = []
for i in range(16):
    significant_QTLs_male = []
    significant_QTLs_female = []
    significant_vQTLs_male = []
    significant_vQTLs_female = []
    GxE_rsIDs_male = []
    GxE_rsIDs_female = []
    for c in chromosomes:

#----------------------------------------------------------------------------------------------------------
# filters all significant QTLs
#----------------------------------------------------------------------------------------------------------

        QTL_path_male = base + "QTL_output/QTL_effects_chr" + c + "_P" + str(i) + "_male.txt"
        QTL_path_female = base + "QTL_output/QTL_effects_chr" + c + "_P" + str(i) + "_female.txt"
        QTLs_male = pd.read_csv(QTL_path_male, delimiter = "\t", header = 0)
        if c != "Y":
            QTLs_female = pd.read_csv(QTL_path_female, delimiter = "\t", header = 0)
        else:
            QTLs_female = COPY(QTLs_male.loc[:, :])
            QTLs_female.loc[:, :] = np.nan

        QTLs_male["chr"] = c
        QTLs_female["chr"] = c

        p_vals_male = QTLs_male["p"].to_numpy()
        p_vals_female = QTLs_female["p"].to_numpy()
        nan_vals_male = np.isnan(p_vals_male)
        nan_vals_male_only = np.logical_and(np.isnan(p_vals_female) == False, np.isnan(p_vals_male))
        nan_vals_female_only = np.logical_and(np.isnan(p_vals_male) == False, np.isnan(p_vals_female))
        chi2_vals = -2*(np.log(p_vals_male) + np.log(p_vals_female))
        p_vals_joint = 1 - chi2.cdf(chi2_vals, df = 4)
        p_vals_joint[nan_vals_male_only] = p_vals_female[nan_vals_male_only]
        p_vals_joint[nan_vals_female_only] = p_vals_male[nan_vals_female_only]

        QTL_is_significant = p_vals_joint < pg
        sig_QTLs_male = QTLs_male.loc[QTL_is_significant, :]
        sig_QTLs_female = QTLs_female.loc[QTL_is_significant, :]
        significant_QTLs_male.append(sig_QTLs_male)
        significant_QTLs_female.append(sig_QTLs_female)
        all_significant_rsIDs_male.append(sig_QTLs_male["rsID"].to_numpy())
        all_significant_rsIDs_female.append(sig_QTLs_female["rsID"].to_numpy())

#----------------------------------------------------------------------------------------------------------
# filters all significant vQTLs
#----------------------------------------------------------------------------------------------------------

        vQTL_path_male = base + "vQTL_output/vQTL_effects_chr" + c + "_P" + str(i) + "_male.txt"
        vQTL_path_female = base + "vQTL_output/vQTL_effects_chr" + c + "_P" + str(i) + "_female.txt"
        vQTLs_male = pd.read_csv(vQTL_path_male, delimiter = "\t", header = 0)
        if c != "Y":
            vQTLs_female = pd.read_csv(vQTL_path_female, delimiter = "\t", header = 0)
        else:
            vQTLs_female = COPY(vQTLs_male.loc[:, :])
            vQTLs_female.loc[:, :] = np.nan

        vQTLs_male["chr"] = c
        vQTLs_female["chr"] = c

        p_vals_male = vQTLs_male["p"].to_numpy()
        p_vals_female = vQTLs_female["p"].to_numpy()
        nan_vals_male_only = np.logical_and(np.isnan(p_vals_female) == False, np.isnan(p_vals_male))
        nan_vals_female_only = np.logical_and(np.isnan(p_vals_male) == False, np.isnan(p_vals_female))
        chi2_vals = -2*(np.log(p_vals_male) + np.log(p_vals_female))
        p_vals_joint = 1 - chi2.cdf(chi2_vals, df = 4)
        p_vals_joint[nan_vals_male_only] = p_vals_female[nan_vals_male_only]
        p_vals_joint[nan_vals_female_only] = p_vals_male[nan_vals_female_only]

        vQTL_is_significant = p_vals_joint < pg
        sig_vQTLs_male = vQTLs_male.loc[vQTL_is_significant, :]
        sig_vQTLs_female = vQTLs_female.loc[vQTL_is_significant, :]
        significant_vQTLs_male.append(sig_vQTLs_male)
        significant_vQTLs_female.append(sig_vQTLs_female)
        all_significant_rsIDs_male.append(sig_vQTLs_male["rsID"].to_numpy())
        all_significant_rsIDs_female.append(sig_vQTLs_female["rsID"].to_numpy())

#----------------------------------------------------------------------------------------------------------
# filters all GxE candidates (rsIDs that are either significant QTLs or vQTLs)
#----------------------------------------------------------------------------------------------------------

        rsID_is_GxE_candidate = np.logical_or(QTL_is_significant, vQTL_is_significant)
        GxE_candidate_rsIDs_male = QTLs_male.loc[rsID_is_GxE_candidate, :]
        GxE_candidate_rsIDs_female = QTLs_female.loc[rsID_is_GxE_candidate, :]
        GxE_rsIDs_male.append(GxE_candidate_rsIDs_male)
        GxE_rsIDs_female.append(GxE_candidate_rsIDs_female)
        all_significant_rsIDs_male.append(sig_QTLs_male["rsID"].to_numpy())
        all_significant_rsIDs_female.append(sig_QTLs_female["rsID"].to_numpy())

    significant_QTLs_male = pd.concat(significant_QTLs_male)
    sig_QTL_path_male = "significant_QTLs/significant_QTLs_P" + str(i) + "_male.txt"
    significant_QTLs_male.to_csv(sig_QTL_path_male, sep = "\t", header = True, index = False)

    significant_QTLs_female = pd.concat(significant_QTLs_female)
    sig_QTL_path_female = "significant_QTLs/significant_QTLs_P" + str(i) + "_female.txt"
    significant_QTLs_female.to_csv(sig_QTL_path_female, sep = "\t", header = True, index = False)

    significant_vQTLs_male = pd.concat(significant_vQTLs_male)
    sig_vQTL_path_male = "significant_vQTLs/significant_vQTLs_P" + str(i) + "_male.txt"
    significant_vQTLs_male.to_csv(sig_vQTL_path_male, sep = "\t", header = True, index = False)

    significant_vQTLs_female = pd.concat(significant_vQTLs_female)
    sig_vQTL_path_female = "significant_vQTLs/significant_vQTLs_P" + str(i) + "_female.txt"
    significant_vQTLs_female.to_csv(sig_vQTL_path_female, sep = "\t", header = True, index = False)

    GxE_rsIDs_male = pd.concat(GxE_rsIDs_male)
    GxE_path_male = "GxE_candidate_rsIDs/GxE_candidate_rsIDs_P" + str(i) + "_male.txt"
    GxE_rsIDs_male.to_csv(GxE_path_male, sep = "\t", header = True, index = False)
    
    GxE_rsIDs_female = pd.concat(GxE_rsIDs_female)
    GxE_path_female = "GxE_candidate_rsIDs/GxE_candidate_rsIDs_P" + str(i) + "_female.txt"
    GxE_rsIDs_female.to_csv(GxE_path_female, sep = "\t", header = True, index = False)

all_significant_rsIDs = pd.DataFrame(reduce(np.union1d, all_significant_rsIDs_male + all_significant_rsIDs_female))
all_significant_rsIDs.to_csv("significant_rsIDs.txt", sep = "\t", header = False, index = False)


