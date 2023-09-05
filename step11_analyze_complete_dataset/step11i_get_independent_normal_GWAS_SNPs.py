import numpy as np
import pandas as pd
import argparse
import os
from functools import reduce
from tqdm import tqdm
from copy import deepcopy as COPY
from scipy.stats import linregress
import statsmodels.api as sm
from scipy.stats import chi2
from bed_reader import open_bed
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
import pdb
import matplotlib as mpl
from scipy.stats import norm
mpl.rcParams['agg.path.chunksize'] = 10000
from statsmodels.stats.outliers_influence import summary_table

base = "/home/greggj/pleiotropy_and_GxE/step9_regress_phenotypes_against_SNPs_logistic_PCA/binary_HF_QTL_output/"
chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
file_name_sets = [[base + name for name in os.listdir(base) if "chr" + i + "_" in name] for i in chromosomes]
files = [pd.concat([pd.read_csv(name, delimiter = "\t") for name in set]) for set in file_name_sets]
for i in range(len(chromosomes)): files[i]["chr"] = chromosomes[i]
files = pd.concat(files)

original_hits = np.array(["rs1906609", "rs7857118", "rs12627426", "rs73839819", "rs2234962", "rs12138073"])
cond1 = files["p_main"] < 5E-8
cond2 = np.isin(files["rsID"].to_numpy().reshape(-1), original_hits)
files = files.loc[np.logical_or(cond1, cond2), :].sort_values(by = "p_main", ascending = True)

def get_conditional_significance(genotypes_raw, offset, phenotype, rsIDs, top_rsID, best_hits, high_LD_hits):

    genotypes = COPY(genotypes_raw).T
    for g in genotypes: g[np.isnan(g)] = np.nanmean(g)
    genotypes = genotypes.T
    # no phenotype imputation because no ICD values are missing (technically, absent = 0 and present = 1)

    k = len(top_rsID)
    pb = 0.05/(len(rsIDs))
    best_hits.append(top_rsID)
    top_indices = np.isin(rsIDs, best_hits)
    top_genotypes = genotypes[:, top_indices]
       
    g_top = genotypes[:, rsIDs == top_rsID].reshape(-1)
    corrs = np.array([nanlinregress(g_top, gi)[2] for gi in genotypes.T])
    high_LD_hits = np.concatenate([high_LD_hits, rsIDs[corrs**2 > 0.8]])
    exited_indices = np.logical_or(top_indices, np.isin(rsIDs, high_LD_hits))

    p_values = []
    for i in range(len(genotypes[0])):
        if exited_indices[i] == True:
            p_values.append(1)
        else:
            g = genotypes[:, i]
            G = np.concatenate([g.reshape(-1,1), top_genotypes], axis = 1)
            X = np.concatenate([G, np.ones((len(G), 1))], axis = 1)
            Z = COPY(X.T)
            for z in Z: z[np.isnan(z)] = np.nanmean(z)
            maxr = np.nanmax(np.corrcoef(Z[:-2])**2 - np.eye(len(Z[:-2])))
            if maxr > 0.8:
                pdb.set_trace()
            X_sub = np.concatenate([top_genotypes, np.ones((len(G), 1))], axis = 1)

            family = sm.genmod.families.Binomial()
            family.link = sm.genmod.families.links.logit()
            model = sm.GLM(phenotype, X, family = family, offset = offset)
            model_sub = sm.GLM(phenotype, X_sub, family = family, offset = offset)
            model_results = model.fit()
            model_sub_results = model_sub.fit()
            LR_test_statistic = 2*(model_results.llf - model_sub_results.llf)
            p_values.append(chi2.sf(LR_test_statistic, 1))
    p_values = np.array(p_values)
    if np.any(p_values <= pb):
        kept_indices = np.logical_or((p_values <= pb), top_indices)
        next_genotypes = genotypes[:, kept_indices]
        next_rsIDs = np.array(rsIDs)[kept_indices]
        next_top_rsID = rsIDs[np.argmin(p_values)]
        
        next_args = [next_genotypes, offset, phenotype, next_rsIDs, next_top_rsID, best_hits, high_LD_hits]
        return(get_conditional_significance(*next_args))
    else:
        return(best_hits)

def nanlinregress(x, y):
    xval_indices = (np.isnan(x) == False)
    yval_indices = (np.isnan(y) == False)
    val_indices = np.logical_and(xval_indices, yval_indices)
    x, y = x[val_indices], y[val_indices]
    return(linregress(x, y))

for chr in np.unique(files["chr"]):

    QTL_eids_path = "../step9_regress_phenotypes_against_SNPs_logistic_PCA/QTL_phenotypes_eids.txt"
    QTL_eids = pd.read_csv(QTL_eids_path, delimiter = "\t", header = None)
    QTL_eids.columns = ["eid"]

    new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
    new_fam_paths = np.array([new_fam_path_prefix + i + ".fam" for i in chromosomes])
    old_fam_path = "../step4_remove_relatives/UKB_samples_unrelated.fam"
    new_fams = [pd.read_csv(path, delim_whitespace = True, header = None)[1].to_numpy() for path in new_fam_paths]
    # new_fam_intersect = reduce(np.intersect1d, new_fams)
    # all eids in QTL_eids exist in new_fam_intersec, but the converse is not true due to modifications. 
    # the new new_fam_intersect term will serve the same function as previously. 
    new_fam = QTL_eids["eid"].to_numpy()
    new_fam_df = pd.DataFrame(np.array([new_fam, np.arange(len(new_fam))]).T)
    new_fam_df.columns = ["eid", "index"] 

    old_fam = pd.read_csv(old_fam_path, delim_whitespace = True, header = None)[0].to_numpy()
    old_in_new = np.isin(old_fam, new_fam)
    if not np.all(old_fam == np.sort(old_fam)):
        print("error1: code expects the fam file from step 4 to have sorted eids")
        exit()
    if not np.all(new_fam == old_fam[old_in_new]):
        print("error2: code expects the fam file from step 4 to have sorted eids")
        exit()

    geno_path = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr" + chr
    pheno_path = "../step7_adjust_HF_for_covariates_logistic_PCA/binary_HF_values.txt"
    PCs_path = "../step7_adjust_HF_for_covariates_logistic_PCA/binary_HF_genetic_PCs.txt"
    Betas_path = "../step7_adjust_HF_for_covariates_logistic_PCA/binary_HF_PC_betas.txt"
    rsIDs = pd.read_csv(geno_path  + ".bim", delim_whitespace = True, header = None)[1]

    phenotypes_info = pd.read_csv(pheno_path, delimiter = "\t").loc[old_in_new, :]
    PCs_info = pd.read_csv(PCs_path, delimiter = "\t").loc[old_in_new, :]
    Betas_info = pd.read_csv(Betas_path, delimiter = "\t")

    prefix9 = "../step9_regress_phenotypes_against_SNPs_logistic_PCA/"
    phenotypes_info_prev = pd.read_csv(prefix9 + "QTL_phenotypes_eids.txt", delimiter = "\t", header = None)
    if not np.all(phenotypes_info["eid"].to_numpy() == phenotypes_info_prev[0].to_numpy()):
        print("exiting: eids are not in the expected order")
        exit()
    if not np.all(PCs_info["eid"].to_numpy() == phenotypes_info_prev[0].to_numpy()):
        print("exiting: eids are not in the expected order")
        exit()

    phenotypes = phenotypes_info.loc[:, "any_HF"].to_numpy()
    Betas = Betas_info.loc[:, "any_HF"].to_numpy()
    PCs = PCs_info[PCs_info.columns[PCs_info.columns != "eid"]].to_numpy()
    PCs2 = np.concatenate([PCs, np.ones((len(PCs), 1))], axis = 1)
    LR_offset_all = np.matmul(PCs2, Betas)

    best_rsIDs = []
 
    main_path = new_fam_paths[np.array(chromosomes) == chr][0]
    new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
    is_in_intersect = np.isin(new_fam_main, new_fam)
    sorted_main_indices = np.argsort(new_fam_main)
    sorted_indices = sorted_main_indices[is_in_intersect[sorted_main_indices]]

    QTLs = files
    QTLs_chr = QTLs.loc[QTLs["chr"].astype(str) == chr, :]

    if "rs1906609" in QTLs_chr["rsID"].to_numpy():
        top_rsID = "rs1906609"
    elif "rs7857118" in QTLs_chr["rsID"].to_numpy():
        top_rsID = "rs7857118"
    elif "rs12627426" in QTLs_chr["rsID"].to_numpy():
        top_rsID = "rs12627426"
    elif "rs73839819" in QTLs_chr["rsID"].to_numpy():
        top_rsID = "rs73839819"
    elif "rs2234962" in QTLs_chr["rsID"].to_numpy():
        top_rsID = "rs2234962"
    elif "rs12138073" in QTLs_chr["rsID"].to_numpy():
        top_rsID = "rs12138073"
    else:
        top_rsID = QTLs_chr["rsID"].to_numpy()[np.argmin(QTLs_chr["p_main"])]

    rsIDs = pd.read_csv(main_path[:-4]  + ".bim", delim_whitespace = True, header = None)[1]
    rsIDs_pheno_indices =  np.where(np.isin(rsIDs, QTLs_chr["rsID"]))[0]
    rsIDs_pheno = rsIDs[rsIDs_pheno_indices]
    genotypes_getter = open_bed(main_path[:-4]  + ".bed", count_A1 = False, num_threads = 1)
    genotypes = genotypes_getter.read(np.s_[sorted_indices,  rsIDs_pheno_indices])
    is_not_maf = np.nanmean(genotypes, axis = 0) > 1 
    genotypes[:, is_not_maf] = 2 - genotypes[:, is_not_maf]
    
    best_rsIDs_pheno = get_conditional_significance(genotypes, LR_offset_all, phenotypes, rsIDs_pheno.to_numpy(), top_rsID, [], np.array([]))
    best_rsIDs_pheno_loc = np.isin(rsIDs_pheno.to_numpy(), best_rsIDs_pheno)
    L = len(best_rsIDs_pheno)
    best_rsID_geno = genotypes[:, best_rsIDs_pheno_loc].T
    for k in range(L): best_rsID_geno[k, np.isnan(best_rsID_geno[k, :])] = np.nanmean(best_rsID_geno[k, :])
    max_R2 = np.max((np.corrcoef(best_rsID_geno)**2 - np.eye(L)))
    best_rsIDs = np.array([best_rsIDs_pheno, L*["binary_AHF"], L*[max_R2]]).T

    best_rsIDs2 = pd.DataFrame(best_rsIDs)
    best_rsIDs2.columns = ["rsID", "pheno_index", "max_LD"]
    path = "normal_GWAS_hits_binary_AHF/QTL_hits_chr" + chr + ".txt"
    if not os.path.isdir("normal_GWAS_hits_binary_AHF"):
        os.mkdir("normal_GWAS_hits_binary_AHF")
    best_rsIDs2.to_csv(path, sep = "\t", header = True, index = False)
