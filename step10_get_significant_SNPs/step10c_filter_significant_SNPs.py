import numpy as np
import pandas as pd
import os
from functools import reduce
from tqdm import tqdm
from scipy.stats import linregress
from bed_reader import open_bed
import pdb

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
pheno_indices = list(range(16))

QTL_pheno_path = "../step9_regress_phenotypes_against_SNPs/QTL_phenotypes.txt"
QTL_phenotypes = pd.read_csv(QTL_pheno_path, delimiter = "\t", header = None).to_numpy()

vQTL_pheno_path = "../step9_regress_phenotypes_against_SNPs/vQTL_phenotypes.txt"
vQTL_phenotypes = pd.read_csv(vQTL_pheno_path, delimiter = "\t", header = None).to_numpy()

QTL_paths = ["significant_QTLs/significant_QTLs_P" + str(i) + ".txt" for i in pheno_indices]
QTLs_info = [pd.read_csv(path, delimiter = "\t") for path in QTL_paths]

vQTL_paths = ["significant_vQTLs/significant_vQTLs_P" + str(i) + ".txt" for i in pheno_indices]
vQTLs_info = [pd.read_csv(path, delimiter = "\t") for path in vQTL_paths]

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

def nanlinregress(x, y):
    xval_indices = (np.isnan(x) == False)
    yval_indices = (np.isnan(y) == False)
    val_indices = np.logical_and(xval_indices, yval_indices)
    x, y = x[val_indices], y[val_indices]
    return(linregress(x, y))

def get_conditional_significance(genotypes, phenotype, rsIDs, top_rsID, best_hits, num_recursions = 0):
    best_hits.append(top_rsID)
    num_recursions += 1
    top_index = np.where(rsIDs == top_rsID)[0][0]
    top_genotype = genotypes[:, top_index]
    slope, intercept = nanlinregress(top_genotype, phenotype)[0:2]
    residuals = phenotype - (slope*top_genotype + intercept)
    new_p_vals = np.array([nanlinregress(genotype, residuals)[3] for genotype in genotypes.T])
    if len(rsIDs) == num_recursions:
        return(best_hits)
    pb = 0.05/(len(rsIDs) - num_recursions)
    if np.any(new_p_vals <= pb):
        # TODO: DOUBLE CHECK THIS WORKS AS EXPECTED
        kept_indices = (new_p_vals <= pb)
        next_genotypes = genotypes[:, kept_indices]
        next_phenotype = residuals
        next_rsIDs = np.array(rsIDs)[kept_indices]
        next_top_rsID = rsIDs[np.argmin(new_p_vals)]
        next_args = [next_genotypes, next_phenotype, next_rsIDs, next_top_rsID, best_hits, num_recursions]
        return(get_conditional_significance(*next_args))
    else:
        return(best_hits)

new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
new_fam_paths = [new_fam_path_prefix + i + ".fam" for i in chromosomes]
old_fam_path = "../step4_remove_relatives/UKB_samples_unrelated.fam"
new_fams = [pd.read_csv(path, delim_whitespace = True, header = None)[1].to_numpy() for path in new_fam_paths]
new_fam_intersect = reduce(np.intersect1d, new_fams)

best_QTLs = []
best_vQTLs = []
for i in pheno_indices:
    QTL_P = QTL_phenotypes[:, i]
    vQTL_P = vQTL_phenotypes[:, i]
    QTLs = QTLs_info[i]
    vQTLs = vQTLs_info[i]
    best_QTLs_pheno = []
    best_vQTLs_pheno = []
    for j in range(len(chromosomes)):
        print(j)
        chr = chromosomes[j]
        path = chr_plink_paths[j]
        
        main_path = new_fam_paths[j]
        new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
        is_in_intersect = np.isin(new_fam_main, new_fam_intersect)
        sorted_main_indices = np.argsort(new_fam_main)
        sorted_indices = sorted_main_indices[is_in_intersect]

        rsIDs = pd.read_csv(path  + ".bim", delim_whitespace = True, header = None)[1]
        genotypes = open_bed(path  + ".bed", count_A1 = False, num_threads = 1).read(np.s_[sorted_indices, :])
        QTLs_chr = QTLs.loc[QTLs["chr"].astype(str) == chr, :]
        vQTLs_chr = vQTLs.loc[vQTLs["chr"].astype(str) == chr, :]
        QTL_rsIDs = QTLs_chr.loc[:, "rsID"].to_numpy()
        vQTL_rsIDs = vQTLs_chr.loc[:, "rsID"].to_numpy()
        QTL_genotypes = genotypes[:, np.isin(rsIDs, QTL_rsIDs)]
        vQTL_genotypes = genotypes[:, np.isin(rsIDs, vQTL_rsIDs)]
        if len(QTL_rsIDs) > 0:
            top_QTL_rsID = QTLs_chr.loc[QTLs_chr["p"] == np.min(QTLs_chr["p"]), "rsID"].values[0]
            best_QTLs_pheno += get_conditional_significance(QTL_genotypes, QTL_P, QTL_rsIDs, top_QTL_rsID, [])
        if len(vQTL_rsIDs) > 0:
            top_vQTL_rsID = vQTLs_chr.loc[vQTLs_chr["p"] == np.min(vQTLs_chr["p"]), "rsID"].values[0]
            best_vQTLs_pheno += get_conditional_significance(vQTL_genotypes, vQTL_P, vQTL_rsIDs, top_vQTL_rsID, [])
    best_QTLs.append(best_QTLs_pheno)
    best_vQTLs.append(best_vQTLs_pheno)

main_QTLs = []
main_vQTLs = []
for i in pheno_indices:
    QTLs = pd.read_csv("significant_QTLs/significant_QTLs_P" + str(i) + ".txt", delimiter = "\t")
    vQTLs = pd.read_csv("significant_vQTLs/significant_vQTLs_P" + str(i) + ".txt", delimiter = "\t")
    QTL_rsids = np.array(best_QTLs[i])
    vQTL_rsids = np.array(best_vQTLs[i])
    QTL_keepers = QTLs.loc[np.isin(QTLs["rsID"], QTL_rsids), :]
    vQTL_keepers = vQTLs.loc[np.isin(vQTLs["rsID"], vQTL_rsids), :]
    del QTL_keepers["slope"]
    del vQTL_keepers["slope"]
    QTL_keepers.columns = ["rsID", "phenotype", "r", "p", "err", "chr"]
    vQTL_keepers.columns = ["rsID", "phenotype", "r", "p", "err", "chr"]
    QTL_keepers.loc[:, "phenotype"] = i
    vQTL_keepers.loc[:, "phenotype"] = i
    main_QTLs.append(QTL_keepers)
    main_vQTLs.append(vQTL_keepers)

pd.concat(main_QTLs).to_csv("QTL_hits.txt", sep = "\t", header = True, index = False)
pd.concat(main_vQTLs).to_csv("vQTL_hits.txt", sep = "\t", header = True, index = False)

