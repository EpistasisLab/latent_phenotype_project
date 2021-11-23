import numpy as np
import pandas as pd
import os
from functools import reduce
from tqdm import tqdm
from copy import deepcopy as COPY
from scipy.stats import linregress
from scipy.stats import chi2
from bed_reader import open_bed
import pdb

if not os.path.exists("GxE_residuals"):
    os.mkdir("GxE_residuals")

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
all_chromosomes = COPY(chromosomes)
pheno_indices = list(range(16))

is_male_path = "../step9_regress_phenotypes_against_SNPs/is_male.txt"
is_male = pd.read_csv(is_male_path, delimiter = "\t", header = None).to_numpy().reshape(-1)

QTL_pheno_path = "../step9_regress_phenotypes_against_SNPs/QTL_phenotypes.txt"
QTL_phenotypes = pd.read_csv(QTL_pheno_path, delimiter = "\t", header = None)
QTL_phenotypes_male = QTL_phenotypes.loc[is_male, :].to_numpy()
QTL_phenotypes_female = QTL_phenotypes.loc[is_male == False, :].to_numpy()

vQTL_pheno_path = "../step9_regress_phenotypes_against_SNPs/vQTL_phenotypes.txt"
vQTL_phenotypes = pd.read_csv(vQTL_pheno_path, delimiter = "\t", header = None)
vQTL_phenotypes_male = vQTL_phenotypes.loc[is_male, :].to_numpy()
vQTL_phenotypes_female = vQTL_phenotypes.loc[is_male == False, :].to_numpy()

env_factors_path = "../step9_regress_phenotypes_against_SNPs/env_factors.txt"
env_factors = pd.read_csv(env_factors_path, delimiter = "\t", header = None)
env_factors_male = env_factors.loc[is_male, :].to_numpy()
env_factors_female = env_factors.loc[is_male == False, :].to_numpy()

QTL_paths_male = ["significant_QTLs/significant_QTLs_P" + str(i) + "_male.txt" for i in pheno_indices]
QTL_paths_female = ["significant_QTLs/significant_QTLs_P" + str(i) + "_female.txt" for i in pheno_indices]
QTLs_info_male = [pd.read_csv(path, delimiter = "\t") for path in QTL_paths_male]
QTLs_info_female = [pd.read_csv(path, delimiter = "\t") for path in QTL_paths_female]

vQTL_paths_male = ["significant_vQTLs/significant_vQTLs_P" + str(i) + "_male.txt" for i in pheno_indices]
vQTL_paths_female = ["significant_vQTLs/significant_vQTLs_P" + str(i) + "_female.txt" for i in pheno_indices]
vQTLs_info_male = [pd.read_csv(path, delimiter = "\t") for path in vQTL_paths_male]
vQTLs_info_female = [pd.read_csv(path, delimiter = "\t") for path in vQTL_paths_female]

GxE_paths_male = ["GxE_candidate_rsIDs/GxE_candidate_rsIDs_P" + str(i) + "_male.txt" for i in pheno_indices]
GxE_paths_female = ["GxE_candidate_rsIDs/GxE_candidate_rsIDs_P" + str(i) + "_female.txt" for i in pheno_indices]
GxE_info_male = [pd.read_csv(path, delimiter = "\t") for path in GxE_paths_male]
GxE_info_female = [pd.read_csv(path, delimiter = "\t") for path in GxE_paths_female]

chr_plink_paths = ["significant_SNPs_plink_files/significant_SNPs_chr" + chr for chr in chromosomes]  
chr_has_hits = []
for path in chr_plink_paths:
    if os.path.isfile(path + ".bed"):
        chr_has_hits.append(True)
    else:
        chr_has_hits.append(False)
chr_has_hits = np.array(chr_has_hits)
chromosomes = np.array(chromosomes)[chr_has_hits]
chr_plink_paths = np.array(chr_plink_paths)[chr_has_hits]

new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
new_fam_paths = np.array([new_fam_path_prefix + i + ".fam" for i in all_chromosomes])
old_fam_path = "../step4_remove_relatives/UKB_samples_unrelated.fam"
new_fams = [pd.read_csv(path, delim_whitespace = True, header = None)[1].to_numpy() for path in new_fam_paths]
new_fam_intersect = reduce(np.intersect1d, new_fams)
new_fam_paths = new_fam_paths[chr_has_hits]

def nanlinregress(x, y):
    xval_indices = (np.isnan(x) == False)
    yval_indices = (np.isnan(y) == False)
    val_indices = np.logical_and(xval_indices, yval_indices)
    x, y = x[val_indices], y[val_indices]
    return(linregress(x, y))

def get_gxE_p_values(genotype, env_factors, pheno_residuals, pb_gxE):

    Y, X = pheno_residuals, np.concatenate((env_factors, np.ones((len(env_factors), 1))), axis = 1)
    XT_X = np.matmul(X.T, X)
    correction = np.eye(len(XT_X))*np.min(XT_X)*1E-6
    XT_X_inv = np.linalg.inv(XT_X + correction)
    coef = np.matmul(XT_X_inv, X.T)
    weights = np.matmul(coef[:, np.isnan(Y) == False], Y[np.isnan(Y) == False])
    Y_est = np.matmul(X, weights)
    residuals2 = Y - Y_est
    p_values = []
    residuals3 = []
    for E in env_factors.T:
        m, b, void, p, void = nanlinregress(genotype*E, residuals2)
        p_values.append(p)
        residuals2_est = m*genotype*E + b
        residuals3.append(residuals2 - residuals2_est)
    return(p_values, residuals3)

def get_GxE_conditional_significance(all_genotypes, QTL_phenotypes, rsIDs,
                                     env_factors, best_hits):

    pb = 0.05/(32*(len(rsIDs)))
    genotypes_male, genotypes_female = all_genotypes
    phenotype_male, phenotype_female = QTL_phenotypes
    env_factors_male, env_factors_female = env_factors
    residuals_male = np.zeros(genotypes_male.shape)
    residuals_female = np.zeros(genotypes_female.shape)
    gxE_pval_arrays_male = []
    gxE_pval_arrays_female = []
    joint_pvals = np.zeros(len(rsIDs))
    for i in range(len(rsIDs)):

        g_male, g_female = COPY(genotypes_male[:, i]), COPY(genotypes_female[:, i])
        slope_male, intercept_male = nanlinregress(g_male, phenotype_male)[0:2]
        slope_female, intercept_female = nanlinregress(g_female, phenotype_female)[0:2]
        res1_male = phenotype_male - (slope_male*g_male + intercept_male)
        res1_female = phenotype_female - (slope_female*g_female + intercept_female)
        gxE_pvals_male, final_residuals_male = get_gxE_p_values(g_male, env_factors_male, res1_male, pb)
        gxE_pvals_female, final_residuals_female = get_gxE_p_values(g_female, env_factors_female, res1_female, pb)
        gxE_pval_arrays_male.append(gxE_pvals_male)
        gxE_pval_arrays_female.append(gxE_pvals_female)
        gender_chi2_vals = -2*(np.log(gxE_pvals_male) + np.log(gxE_pvals_female))
        gxE_pvals = [chi2.sf(val, df = 4) for val in gender_chi2_vals]
        joint_pvals[i] = np.min(gxE_pvals)
        residuals_male[:, i] += final_residuals_male[np.argmin(gxE_pvals)]
        residuals_female[:, i] += final_residuals_female[np.argmin(gxE_pvals)]

    print("pb: " + str(pb))
    print("min p: " + str(np.min(joint_pvals)))
    if np.any(joint_pvals <= pb):
        print(joint_pvals[joint_pvals <= pb])
        ind = np.where(joint_pvals == np.min(joint_pvals))[0][0]
        data = [rsIDs[ind], joint_pvals[ind]]
        data += gxE_pval_arrays_male[ind] + gxE_pval_arrays_female[ind]
        best_hits.append(data)
        kept_indices = (joint_pvals <= pb)
        next_genotypes = [genotypes_male[:, kept_indices], genotypes_female[:, kept_indices]]
        next_phenotypes = [residuals_male[:, ind], residuals_female[:, ind]]
        env_factors = [env_factors_male, env_factors_female]
        next_rsIDs = rsIDs[kept_indices]
        next_args = [next_genotypes, next_phenotypes, next_rsIDs]
        next_args += [env_factors, best_hits]
        return(get_GxE_conditional_significance(*next_args))
    else:
        return(best_hits, QTL_phenotypes)

def get_conditional_significance(genotypes_male, genotypes_female, phenotype_male, phenotype_female, 
                                 rsIDs, top_rsID, best_hits):

    best_hits.append(top_rsID)
    pb = 0.05/(len(rsIDs))
    top_index = np.where(rsIDs == top_rsID)[0][0]
    top_genotype_male = genotypes_male[:, top_index]
    top_genotype_female = genotypes_female[:, top_index]
    slope_male, intercept_male = nanlinregress(top_genotype_male, phenotype_male)[0:2]
    slope_female, intercept_female = nanlinregress(top_genotype_female, phenotype_female)[0:2]
    residuals_male = phenotype_male - (slope_male*top_genotype_male + intercept_male)
    residuals_female = phenotype_female - (slope_female*top_genotype_female + intercept_female)
    p_vals_male = np.array([nanlinregress(genotype, residuals_male)[3] for genotype in genotypes_male.T])
    p_vals_female = np.array([nanlinregress(genotype, residuals_female)[3] for genotype in genotypes_female.T])
    nan_vals_male_only = np.logical_and(np.isnan(p_vals_female) == False, np.isnan(p_vals_male))
    nan_vals_female_only = np.logical_and(np.isnan(p_vals_male) == False, np.isnan(p_vals_female))
    chi2_vals = -2*(np.log(p_vals_male) + np.log(p_vals_female))
    p_vals_joint = 1 - chi2.cdf(chi2_vals, df = 4)
    p_vals_joint[nan_vals_male_only] = p_vals_female[nan_vals_male_only]
    p_vals_joint[nan_vals_female_only] = p_vals_male[nan_vals_female_only]
    if np.any(p_vals_joint <= pb):
        kept_indices = (p_vals_joint <= pb)
        next_genotypes_male = genotypes_male[:, kept_indices]
        next_genotypes_female = genotypes_female[:, kept_indices]
        next_phenotype_male = residuals_male
        next_phenotype_female = residuals_female
        next_rsIDs = np.array(rsIDs)[kept_indices]
        next_top_rsID = rsIDs[np.argmin(p_vals_joint)]
        next_args = [next_genotypes_male, next_genotypes_female, next_phenotype_male, next_phenotype_female]
        next_args += [next_rsIDs, next_top_rsID, best_hits]
        return(get_conditional_significance(*next_args))
    else:
        return(best_hits)

best_QTLs = []
best_vQTLs = []
best_GxE = []
male_residuals = []
female_residuals = []
for i in pheno_indices:
    print("i: " + str(i))
    QTL_P_male = QTL_phenotypes_male[:, i]
    QTL_P_female = QTL_phenotypes_female[:, i]
    vQTL_P_male = vQTL_phenotypes_male[:, i]
    vQTL_P_female = vQTL_phenotypes_female[:, i]
    QTLs_male = QTLs_info_male[i]
    QTLs_female = QTLs_info_female[i]
    vQTLs_male = vQTLs_info_male[i]
    vQTLs_female = vQTLs_info_female[i]
    GxE_male = GxE_info_male[i]
    GxE_female = GxE_info_female[i]
    best_QTLs_pheno = []
    best_vQTLs_pheno = []
    best_GxE_pheno = []
    male_residuals_pheno = []
    female_residuals_pheno = []
    for j in range(len(chromosomes)):
        print(j)
        if str(j) == "22":
            pdb.set_trace() 
        chr = chromosomes[j]
        path = chr_plink_paths[j]
        
        main_path = new_fam_paths[j]
        new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
        is_in_intersect = np.isin(new_fam_main, new_fam_intersect)
        sorted_main_indices = np.argsort(new_fam_main)
        sorted_indices = sorted_main_indices[is_in_intersect[sorted_main_indices]]

        rsIDs = pd.read_csv(path  + ".bim", delim_whitespace = True, header = None)[1]
        genotypes = open_bed(path  + ".bed", count_A1 = False, num_threads = 1).read()
        genotypes = genotypes[sorted_indices, :]
        is_not_maf = np.nanmean(genotypes, axis = 0) > 1 
        genotypes[:, is_not_maf] = 2 - genotypes[:, is_not_maf]

        QTLs_male_chr = QTLs_male.loc[QTLs_male["chr"].astype(str) == chr, :]
        QTLs_female_chr = QTLs_female.loc[QTLs_female["chr"].astype(str) == chr, :]
        
        # all of the same rsIDs are examined for both males and females.
        QTL_rsIDs = QTLs_male_chr.loc[:, "rsID"].to_numpy()
        QTL_rsIDs = rsIDs[np.isin(rsIDs, QTL_rsIDs)].to_numpy()
        QTL_genotypes = genotypes[:, np.isin(rsIDs, QTL_rsIDs)]

        QTL_male_genotypes = QTL_genotypes[is_male]
        QTL_female_genotypes = QTL_genotypes[is_male == False]
        QTL_male_alphas =  QTLs_male_chr["alpha"].to_numpy()
        QTL_female_alphas =  QTLs_female_chr["alpha"].to_numpy()
        for k in range(len(QTL_male_alphas)): 
            if QTL_male_alphas[k] != -1:
                QTL_male_genotypes[QTL_male_genotypes[:, k] == 1, k] = QTL_male_alphas[k]
        for k in range(len(QTL_female_alphas)): 
            if QTL_female_alphas[k] != -1:
                QTL_female_genotypes[QTL_female_genotypes[:, k] == 1, k] = QTL_female_alphas[k]
        QTL_male_genotypes[QTL_male_genotypes == 2] = 1
        QTL_female_genotypes[QTL_female_genotypes == 2] = 1

        vQTLs_male_chr = vQTLs_male.loc[vQTLs_male["chr"].astype(str) == chr, :]
        vQTLs_female_chr = vQTLs_female.loc[vQTLs_female["chr"].astype(str) == chr, :]
        
        # all of the same rsIDs are examined for both males and females.
        vQTL_rsIDs = vQTLs_male_chr.loc[:, "rsID"].to_numpy()
        vQTL_rsIDs = rsIDs[np.isin(rsIDs, vQTL_rsIDs)].to_numpy()
        vQTL_genotypes = genotypes[:, np.isin(rsIDs, vQTL_rsIDs)]

        vQTL_male_genotypes = vQTL_genotypes[is_male]
        vQTL_female_genotypes = vQTL_genotypes[is_male == False]
        vQTL_male_alphas =  vQTLs_male_chr["alpha"].to_numpy()
        vQTL_female_alphas =  vQTLs_female_chr["alpha"].to_numpy()

        for k in range(len(vQTL_male_alphas)): vQTL_male_genotypes[vQTL_male_genotypes[:, k] == 1, k] = vQTL_male_alphas[k]
        for k in range(len(vQTL_female_alphas)): vQTL_female_genotypes[vQTL_female_genotypes[:, k] == 1, k] = vQTL_female_alphas[k]
        vQTL_male_genotypes[vQTL_male_genotypes == 2] = 1
        vQTL_female_genotypes[vQTL_female_genotypes == 2] = 1

        GxE_male_chr = GxE_male.loc[GxE_male["chr"].astype(str) == chr, :]
        GxE_female_chr = GxE_female.loc[GxE_female["chr"].astype(str) == chr, :]
        
        # all of the same rsIDs are examined for both males and females.
        GxE_rsIDs = GxE_male_chr.loc[:, "rsID"].to_numpy()
        GxE_rsIDs = rsIDs[np.isin(rsIDs, GxE_rsIDs)].to_numpy()
        GxE_genotypes = genotypes[:, np.isin(rsIDs, GxE_rsIDs)]

        GxE_male_genotypes = GxE_genotypes[is_male]
        GxE_female_genotypes = GxE_genotypes[is_male == False]
        GxE_male_alphas =  GxE_male_chr["alpha"].to_numpy()
        GxE_female_alphas =  GxE_female_chr["alpha"].to_numpy()
        for k in range(len(GxE_male_alphas)): 
            if GxE_male_alphas[k] != -1:
                GxE_male_genotypes[GxE_male_genotypes[:, k] == 1, k] = GxE_male_alphas[k]
        for k in range(len(GxE_female_alphas)): 
            if GxE_female_alphas[k] != -1:
                GxE_female_genotypes[GxE_female_genotypes[:, k] == 1, k] = GxE_female_alphas[k]
        GxE_male_genotypes[GxE_male_genotypes == 2] = 1
        GxE_female_genotypes[GxE_female_genotypes == 2] = 1


        if len(QTL_rsIDs) > 0:

            male_nans = np.isnan(QTL_male_alphas)
            female_nans = np.isnan(QTL_female_alphas)
            QTL_male_genotypes[:, male_nans] = np.nan
            QTL_female_genotypes[:, female_nans] = np.nan

            p_vals_male = QTLs_male_chr["p"].to_numpy()
            p_vals_female = QTLs_female_chr["p"].to_numpy()
            nan_vals_male_only = np.logical_and(np.isnan(p_vals_female) == False, np.isnan(p_vals_male))
            nan_vals_female_only = np.logical_and(np.isnan(p_vals_male) == False, np.isnan(p_vals_female))
            chi2_vals = -2*(np.log(p_vals_male) + np.log(p_vals_female))
            p_vals_joint = 1 - chi2.cdf(chi2_vals, df = 4)
            p_vals_joint[nan_vals_male_only] = p_vals_female[nan_vals_male_only]
            p_vals_joint[nan_vals_female_only] = p_vals_male[nan_vals_female_only]
            top_QTL_rsID = QTLs_male_chr.loc[p_vals_joint == np.min(p_vals_joint), "rsID"].values[0]
            best_QTLs_pheno += get_conditional_significance(QTL_male_genotypes, QTL_female_genotypes, 
                                                            QTL_P_male, QTL_P_female, QTL_rsIDs, top_QTL_rsID, [])
        if len(vQTL_rsIDs) > 0:

            male_nans = np.isnan(vQTL_male_alphas)
            female_nans = np.isnan(vQTL_female_alphas)
            vQTL_male_genotypes[:, male_nans] = np.nan
            vQTL_female_genotypes[:, female_nans] = np.nan

            p_vals_male = vQTLs_male_chr["p"].to_numpy()
            p_vals_female = vQTLs_female_chr["p"].to_numpy()
            nan_vals_male_only = np.logical_and(np.isnan(p_vals_female) == False, np.isnan(p_vals_male))
            nan_vals_female_only = np.logical_and(np.isnan(p_vals_male) == False, np.isnan(p_vals_female))
            chi2_vals = -2*(np.log(p_vals_male) + np.log(p_vals_female))
            p_vals_joint = 1 - chi2.cdf(chi2_vals, df = 4)
            p_vals_joint[nan_vals_male_only] = p_vals_female[nan_vals_male_only]
            p_vals_joint[nan_vals_female_only] = p_vals_male[nan_vals_female_only]
            top_vQTL_rsID = vQTLs_male_chr.loc[p_vals_joint == np.min(p_vals_joint), "rsID"].values[0]
            best_vQTLs_pheno += get_conditional_significance(vQTL_male_genotypes, vQTL_female_genotypes, 
                                                             vQTL_P_male, vQTL_P_female, vQTL_rsIDs, top_vQTL_rsID, [])

        if len(GxE_rsIDs) > 0:
            
            male_nans = np.isnan(GxE_male_alphas)
            female_nans = np.isnan(GxE_female_alphas)
            GxE_male_genotypes[:, male_nans] = np.nan
            GxE_female_genotypes[:, female_nans] = np.nan
            genotypes = [GxE_male_genotypes, GxE_female_genotypes]
            phenotypes = [QTL_P_male, QTL_P_female]
            env_factors = [env_factors_male, env_factors_female]
            best_GxE_pheno_subset, GxE_residuals = get_GxE_conditional_significance(genotypes, phenotypes, GxE_rsIDs, env_factors, [])
            best_GxE_pheno += best_GxE_pheno_subset
            male_residuals_pheno.append(GxE_residuals[0])
            female_residuals_pheno.append(GxE_residuals[1])

    best_QTLs.append(best_QTLs_pheno)
    best_vQTLs.append(best_vQTLs_pheno)
    best_GxE.append(best_GxE_pheno)
    male_residuals.append(np.array(male_residuals_pheno))
    female_residuals.append(np.array(female_residuals_pheno))


for gender in ["_male", "_female"]:
    main_QTLs = []
    main_vQTLs = []
    main_GxE = []
    for i in pheno_indices[:-1]:
        QTLs = pd.read_csv("significant_QTLs/significant_QTLs_P" + str(i) + gender + ".txt", delimiter = "\t")
        vQTLs = pd.read_csv("significant_vQTLs/significant_vQTLs_P" + str(i) + gender + ".txt", delimiter = "\t")
        GxE = pd.read_csv("GxE_candidate_rsIDs/GxE_candidate_rsIDs_P" + str(i) + gender + ".txt", delimiter = "\t")

        QTL_rsids = np.array(best_QTLs[i])
        vQTL_rsids = np.array(best_vQTLs[i])
        GxE_rsids = pd.DataFrame(best_GxE[i])
        GxE_rsids_len = len(GxE_rsids)
        if GxE_rsids_len == 0:
            GxE_rsids = pd.DataFrame([68*[np.nan]])
        colnames = ["rsID", "jointp"] 
        colnames += ["GxE" + str(i) + "_p_m" for i in range(1, 34)] 
        colnames += ["GxE" + str(i) + "_p_f" for i in range(1, 34)]
        GxE_rsids.columns = colnames
        if GxE_rsids_len == 0: 
            GxE_rsids = GxE_rsids.dropna()

        QTL_keepers = QTLs.loc[np.isin(QTLs["rsID"], QTL_rsids), :]
        vQTL_keepers = vQTLs.loc[np.isin(vQTLs["rsID"], vQTL_rsids), :]
        GxE_keepers = GxE.merge(GxE_rsids, on = "rsID", how = "inner")

        del QTL_keepers["slope"]
        del vQTL_keepers["slope"]
        del GxE_keepers["slope"]

        QTL_keepers.columns = ["rsID", "phenotype", "r", "p", "err", "alpha", "chr"]
        vQTL_keepers.columns = ["rsID", "phenotype", "r", "p", "err", "alpha", "chr"]
        GxE_keepers.columns = ["rsID", "phenotype", "r", "p", "err", "alpha", "chr"] + colnames[1:]

        QTL_keepers.loc[:, "phenotype"] = i
        vQTL_keepers.loc[:, "phenotype"] = i
        GxE_keepers.loc[:, "phenotype"] = i

        main_QTLs.append(QTL_keepers)
        main_vQTLs.append(vQTL_keepers)
        main_GxE.append(GxE_keepers)

    pd.concat(main_QTLs).to_csv("QTL_hits" + gender + ".txt", sep = "\t", header = True, index = False)
    pd.concat(main_vQTLs).to_csv("vQTL_hits" + gender + ".txt", sep = "\t", header = True, index = False)
    pd.concat(main_GxE).to_csv("GxE_hits" + gender + ".txt", sep = "\t", header = True, index = False)
