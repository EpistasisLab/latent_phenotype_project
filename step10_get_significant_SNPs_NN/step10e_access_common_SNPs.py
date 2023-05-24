import numpy as np
import pandas as pd
import statsmodels.api as sm
from itertools import combinations
from functools import reduce
from copy import deepcopy as COPY
from bed_reader import open_bed
from matplotlib import pyplot as plt
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import wilcoxon
from scipy.stats import rankdata
from scipy.stats import mannwhitneyu as mwu
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier as GBC
from sklearn.model_selection import KFold
from scipy.stats import chi2
from scipy.stats import combine_pvalues
from scipy.stats import cauchy
import statsmodels.api as sm
from tqdm import tqdm
import os
import pdb

def cross_validate(X, y, k, model, model2 = None):
    kf = KFold(n_splits = k, shuffle = True)
    kf.get_n_splits(X)
    
    if model2 == None:
        r2 = []
        for train_ind, test_ind in kf.split(X):
            model_info = (COPY(model)).fit(X[train_ind], y[train_ind])
            pps = model_info.predict_proba(X[test_ind])[:, 1]
            r2.append(pearsonr(pps, y[test_ind])[0]**2)
        return(r2)
    else:
        r21_set, r22_set = [], []
        for train_ind, test_ind in kf.split(X):
            model_info1 = (COPY(model)).fit(X[train_ind], y[train_ind])
            pps1 = model_info1.predict_proba(X[test_ind])[:, 1]
            r21_set.append(pearsonr(pps1, y[test_ind])[0]**2)

            model_info2 = (COPY(model2)).fit(X[train_ind], y[train_ind])
            pps2 = model_info2.predict_proba(X[test_ind])[:, 1]
            r22_set.append(pearsonr(pps2, y[test_ind])[0]**2)

        return(np.array(r21_set), np.array(r22_set))



def transform_genotypes(genotypes, phenotypes, Env_factor, impute = True):
    
    g_new, p, E1 = np.zeros(genotypes.shape), COPY(phenotypes), COPY(Env_factor)
        
    for i in range(len(genotypes[0])):

        g = COPY(genotypes[:, i])
        if impute:
            g[np.isnan(g)] = np.nanmedian(g)
        else:
            val_indices = np.logical_or(np.isnan(g), np.isnan(p)) == False
            g, p, E1 = g[val_indices], p[val_indices], E1[val_indices]

        g2_vals = (g == 2)
        if np.all(np.sum(g2_vals, axis = 0) > 1000) and chr != "Y": 

            g1_vals = (g == 1)
            X = np.array([g1_vals, g2_vals, (g1_vals)*E1, (g2_vals)*E1, E1]).T
            X2 = np.array([g1_vals, g2_vals]).T
            y = (p - np.mean(p))

            Betas = LinearRegression().fit(X, y).coef_
            enc1, enc2 = Betas[0] + Betas[2]*E1[g1_vals], Betas[1] + Betas[3]*E1[g2_vals]  
            g_old = COPY(g)
            g[g1_vals] = enc1     
            g[g2_vals] = enc2   
            ge = np.array([g, E1, np.ones(len(g))]).T

        else:
            g[g2_vals] = 1
            g1_vals = (g == 1)
            X = np.array([g1_vals, (g1_vals)*E1, E1]).T
            y = (p - np.mean(p))

            Betas = LinearRegression().fit(X, y).coef_
            enc = Betas[0] + Betas[1]*E1[g1_vals]  
            g_old = COPY(g)
            g[g1_vals] = enc 

        if impute:
            g_new[:, i] = g
        else:
            g_new[val_indices, i] = g
            g_new[val_indices == False, i] = np.nan

    return(g_new)

def nanlinregress(x, y):
    xval_indices = (np.isnan(x) == False)
    yval_indices = (np.isnan(y) == False)
    val_indices = np.logical_and(xval_indices, yval_indices)
    x, y = x[val_indices], y[val_indices]
    return(linregress(x, y))

GxE_effects_fnames = [name for name in os.listdir("hits_GxE_p_vals")]
GxE_effects_info = [name[:-4].split("_") for name in GxE_effects_fnames]
chr = [info[-1] for info in GxE_effects_info]
pheno_index = [int(info[-2]) for info in GxE_effects_info]
env_factor = [info[-3] for info in GxE_effects_info]
rsID = ["_".join(info[:-3]) for info in GxE_effects_info]
GxE_effects_pvals = [pd.read_csv("hits_GxE_p_vals/" + name, sep = "\t", header = None) for name in GxE_effects_fnames]
GxE_effects_pvals = [df.to_numpy()[0].tolist() for df in GxE_effects_pvals]
pvals0 = [i[0] for i in GxE_effects_pvals]
pvals1 = [i[1] for i in GxE_effects_pvals]
pvals2 = [i[2] for i in GxE_effects_pvals]
pvals_complete_set = np.array([i for i in zip(pvals0, pvals1)])
GxE_effects_df = pd.DataFrame(zip(rsID, env_factor, pheno_index, chr, pvals0, pvals1, pvals2))
GxE_effects_df.columns = ["rsID", "env_factor", "pheno_index", "chr", "pEDGE2", "pEDGE", "p_joined"]

# formula for corrected m (g*) in https://pubmed.ncbi.nlm.nih.gov/22588983/
# formula for ICC: https://en.wikipedia.org/wiki/Intraclass_correlation
pheno_path = "../step9_regress_phenotypes_against_SNPs_NN/QTL_phenotypes.txt"
X = pd.read_csv(pheno_path, delimiter = "\t", header = None).to_numpy()[:, :15]
N, K = len(X), len(X[0])
xn, x, s2 = np.mean(X, axis = 1), np.mean(X), np.var(X)
ICC = (K*np.sum((xn - x)**2)/(N*(K - 1)*s2)) - (1/(K - 1))
m = (K + 1) - (1 + (K - 1)*ICC)
pb = 5E-8/(m*6)

# getting main effects with classical testing
suffixes = ["smoking", "alcohol", "exercise", "gender"]
prefixes = ["hits_QTL_" + suf for suf in suffixes]
all_paths = []
for prefix in prefixes: all_paths += [prefix + "/" + path for path in os.listdir(prefix)]
env_factors = [path.split("/")[0].split("_")[-1] for path in all_paths]
all_main_effects = [pd.read_csv(path, delimiter = "\t") for path in all_paths]
for i in range(len(all_main_effects)): all_main_effects[i]["env_factor"] = env_factors[i]
all_main_effects = pd.concat(all_main_effects)
main_effects = (all_main_effects[all_main_effects["p_null2"] <= pb])
main_effects = main_effects.drop_duplicates("rsID")

# getting permutation effects
GxE_effects = GxE_effects_df[GxE_effects_df["pEDGE2"] <= pb].sort_values(by = "pEDGE2")
EDGE_effects = GxE_effects_df[GxE_effects_df["pEDGE"] <= pb].sort_values(by = "pEDGE")
mixed_inds = np.logical_and(GxE_effects_df["p_joined"] <= pb, GxE_effects_df["pEDGE"] > pb)
mixed_inds = np.logical_and(mixed_inds, GxE_effects_df["pEDGE2"] > pb)
mixed_effects = GxE_effects_df[mixed_inds].sort_values(by = "p_joined")

#-----------------------------------------------------------------------------------------------------------------------
# SANITY CHECK
# this part only shows that EDGE and main effect p values are strongly correlated (r = 0.89)
#-----------------------------------------------------------------------------------------------------------------------
main_EDGE_effects = all_main_effects.merge(EDGE_effects, on = ["rsID", "pheno_index", "env_factor"], how = "inner")
p_null2, p_EDGE = np.log10(main_EDGE_effects["p_null2"]), np.log10(main_EDGE_effects["pEDGE"])
non_subbed_inds = (np.abs(p_null2 - p_EDGE) < 1E-6) == False
corr, corrp = spearmanr(p_null2[non_subbed_inds], p_EDGE[non_subbed_inds])
#-----------------------------------------------------------------------------------------------------------------------
# SANITY CHECK END
#-----------------------------------------------------------------------------------------------------------------------

main_effect_counts = [np.sum(main_effects["pheno_index"] == i) for i in range(15)]
GxE_effect_counts = []
column_names = ["phenotype_index", "main effects", "smoking", "alcohol", "gender", "exercise"]
for eff_name in column_names[2:]:
    GxE_eff_name = GxE_effects[GxE_effects["env_factor"] == eff_name]
    GxE_effect_counts.append([np.sum(GxE_eff_name["pheno_index"] == i) for i in range(15)])
step10e_count_data = [[i + 1 for i in range(15)]] + [main_effect_counts] + GxE_effect_counts
step10e_count_data = pd.DataFrame(np.array(step10e_count_data).T)
step10e_count_data.columns = column_names
step10e_count_data.to_csv("step10e_count_data.txt", sep = "\t", header = True, index = False)

Main = main_effects[["rsID", "p_null2", "pheno_index"]]
GxSmoking = GxE_effects.loc[GxE_effects["env_factor"] == "smoking", ["rsID", "pEDGE2", "pheno_index"]]
GxAlcohol = GxE_effects.loc[GxE_effects["env_factor"] == "alcohol", ["rsID", "pEDGE2", "pheno_index"]]
GxGender = GxE_effects.loc[GxE_effects["env_factor"] == "gender", ["rsID", "pEDGE2", "pheno_index"]]
GxExercise = GxE_effects.loc[GxE_effects["env_factor"] == "exercise", ["rsID", "pEDGE2", "pheno_index"]]
Main.to_csv("rsIDs_main_effects.txt", sep = "\t", header = True, index = False)
GxSmoking.to_csv("rsIDs_GxSmoking_effects.txt", sep = "\t", header = True, index = False)
GxAlcohol.to_csv("rsIDs_GxAlcohol_effects.txt", sep = "\t", header = True, index = False)
GxGender.to_csv("rsIDs_GxGender_effects.txt", sep = "\t", header = True, index = False)
GxExercise.to_csv("rsIDs_GxExercise_effects.txt", sep = "\t", header = True, index = False)

eid_path = "../step9_regress_phenotypes_against_SNPs_NN/QTL_phenotypes_eids.txt"
eids = pd.read_csv(eid_path, delimiter = "\t", header = None)
eids.columns = ["eid"]

names = ["_alcohol", "_smoking", "_gender"]
env_IDs = ['annual-consumption', 'pack-years', '22001-0.0']
exercise_IDs = ['874-average', '894-average', '914-average']
path = "../step7_adjust_HF_for_covariates_NN/env_factors_for_step9.txt"
env_data = pd.read_csv(path, delimiter = "\t", usecols = ["eid"] + env_IDs + exercise_IDs)

exercise = env_data[exercise_IDs].to_numpy()
pca = PCA(n_components = 1)
pca.fit(exercise)
exercise_score = pca.transform(exercise).reshape(-1)
env_data["exercise_score"] = exercise_score
env_data = env_data[["eid"] + env_IDs + ["exercise_score"]]
env_data = env_data.merge(eids, on = "eid", how = "inner")
del env_data["eid"]
env_data = env_data.to_numpy()

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
new_fam_paths = np.array([new_fam_path_prefix + i + ".fam" for i in chromosomes])
old_fam_path = "../step4_remove_relatives/UKB_samples_unrelated.fam"
new_fams = [pd.read_csv(path, delim_whitespace = True, header = None)[1].to_numpy() for path in new_fam_paths]
new_fam_intersect = reduce(np.intersect1d, new_fams)

prefixes = ["hits_QTL" + name + "/" for name in names]
path_sets = [[prefix + path for path in os.listdir(prefix)] for prefix in prefixes]
files = [pd.concat([pd.read_csv(path, delimiter = "\t") for path in set]) for set in path_sets]

plink_prefixes = ["significant_SNPs_plink_files" + name + "/" for name in names]
plink_path_sets = [[prefix + path for path in os.listdir(prefix)] for prefix in plink_prefixes]
plink_bed_sets = [[path for path in set if ".bed" in path] for set in plink_path_sets]
plink_bim_sets = [[path for path in set if ".bim" in path] for set in plink_path_sets]
plink_fam_sets = [[path for path in set if ".fam" in path] for set in plink_path_sets]
plink_chr_sets = [np.array([path.split("chr")[1].split(".")[0] for path in set]) for set in plink_bed_sets]

chromosomes_used = pd.concat([file[["chr"]] for file in files])["chr"]
chromosomes_used = np.unique(chromosomes_used.to_numpy(dtype = str))
phenotypes = pd.concat([file[["pheno_index"]] for file in files])["pheno_index"]
phenotypes = np.unique(phenotypes.to_numpy())

if not os.path.isfile("dummy.bim"):
    dummy = pd.DataFrame(2*[6*["dummy"]])
    dummy.to_csv("dummy.bim", sep = "\t", header = False, index = False)
if not os.path.isfile("dummy.bed"):
    file = open("dummy.bed", "w")
    file.close()

col_names = ["chr", "pheno_index", "env1", "env2", "rsid1", "rsid2", "r2"]
col_data = []
EDGE2_geno_pheno_sets = [[] for i in range(len(phenotypes))]
geno_pheno_sets = [[] for i in range(len(phenotypes))]
standard_p_values_rsids = []
standard_p_values_E_ind = []
standard_p_values_GxE = []
figure_geno = 0
figure_geno_EDGE = 0
figure_pheno = 0
for chr in tqdm(chromosomes_used):

    chr_index = np.where(chr == np.array(chromosomes))[0][0]
    main_path = new_fam_paths[chr_index]
    new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
    is_in_intersect = np.isin(new_fam_main, eids.to_numpy().reshape(-1))
    sorted_main_indices = np.argsort(new_fam_main)
    sorted_indices = sorted_main_indices[is_in_intersect[sorted_main_indices]]

    main_effects_chr = main_effects.loc[main_effects["chr"] == chr, ["rsID", "pheno_index"]]
    prefix = "significant_SNPs_plink_files_gender/"
    main_effects_bim_path = prefix + "significant_SNPs_chr" + chr + ".bim"
    if os.path.isfile(main_effects_bim_path):
        main_effects_rsIDs = pd.read_csv(main_effects_bim_path, delim_whitespace = True, header = None)
        main_effects_bed_path = prefix + "significant_SNPs_chr" + chr + ".bed"
        main_effects_bed_opener = open_bed(main_effects_bed_path, count_A1 = False, num_threads = 1)
        for i in phenotypes:
            main_effects_chr_pheno = main_effects_chr.loc[main_effects_chr["pheno_index"] == i, "rsID"]
            main_effects_ind = np.where(main_effects_rsIDs.isin(main_effects_chr_pheno.to_numpy()))[0]
            main_effects_geno = main_effects_bed_opener.read(np.s_[sorted_indices,  main_effects_ind]).T
            for g in main_effects_geno:
                g[np.isnan(g)] = np.nanmedian(g)
                if np.mean(g) > 1: g = 2 - g
                if np.sum(g == 2) < 1000: g[g == 2] = 1
                EDGE2_geno_pheno_sets[i].append(g)
                geno_pheno_sets[i].append(g)
    
    chr_indices = [np.where(set == chr)[0] for set in plink_chr_sets]
    all_E_names, num_E = [name[1:] for name in names], range(len(chr_indices))
    E_name_index_pairs = combinations(num_E, 2)

    QTL_pheno_path = "../step9_regress_phenotypes_against_SNPs_NN/QTL_phenotypes.txt"
    QTL_phenotypes = pd.read_csv(QTL_pheno_path, delimiter = "\t", header = None).to_numpy()
    QTL_phenotypes = (QTL_phenotypes - np.nanmean(QTL_phenotypes, axis = 0))/(np.nanstd(QTL_phenotypes, axis = 0))
    
    #for z in range(16):
    #    p1, p99 = np.percentile(QTL_phenotypes[:, z], [1, 99])
    #    QTL_phenotypes[QTL_phenotypes[:, z] < p1, z] = p1 
    #    QTL_phenotypes[QTL_phenotypes[:, z] > p99, z] = p99

    chr_bed_paths = [np.array(plink_bed_sets[k])[chr_indices[k]] for k in num_E]
    chr_bed_paths = ['' if len(path) == 0 else path[0] for path in chr_bed_paths]
    chr_bim_paths = [np.array(plink_bim_sets[k])[chr_indices[k]] for k in num_E]
    chr_bim_paths = ["dummy.bim" if len(path) == 0 else path[0] for path in chr_bim_paths]

    chr_bed_openers = ['' if len(path) == 0 else open_bed(path, count_A1 = False, num_threads = 1) for path in chr_bed_paths]
    chr_bim_files = [pd.read_csv(path, delim_whitespace = True, header = None) for path in chr_bim_paths]
    results_partial = GxE_effects[GxE_effects["chr"] == chr]   

    for i in phenotypes:
        p = COPY(QTL_phenotypes[:, i])
        p2 = COPY(QTL_phenotypes[:, i])
        p[np.isnan(p)] = np.nanmedian(p)
        geno_pheno_sets 
        results =  results_partial[results_partial["pheno_index"] == i]
        E_names = np.unique(results["env_factor"])
        if len(results) > 0:
            hits_subsets = [results.loc[results["env_factor"] == name, "rsID"] for name in all_E_names]
            hits_indices = [np.where(chr_bim_files[k][1].isin(hits_subsets[k].to_numpy()))[0] for k in num_E]
            genotypes = [np.array([]) if len(hits_indices[k]) == 0 else chr_bed_openers[k].read(np.s_[:,  hits_indices[k]]).T for k in num_E]
            if len(E_names) > 1:
                for i1, i2 in  E_name_index_pairs:
                    g1, g2 = genotypes[i1], genotypes[i2]
                    r1, r2 = range(len(g1)), range(len(g2))
                    for i3 in r1:
                        for i4 in r2:
                            r = nanlinregress(g1[i3], g2[i4])[3]
                            if r**2 > 0.7:
                                E1, E2 = all_E_names[i1], all_E_names[i2]
                                rsID1 = hits_subsets[i1].to_numpy()[i3]
                                rsID2 = hits_subsets[i2].to_numpy()[i4]
                                col_data.append([chr, i, E1, E2, rsID1, rsID2, r**2])

            for j in num_E:
                rsIDs_j = hits_subsets[j].to_numpy()
                geno_sub = genotypes[j]
                if len(geno_sub) > 0:
                    geno_sub =  geno_sub[:, sorted_indices].T
                    afs = np.nanmean(geno_sub, axis = 0)
                    af_is_not_maf = (afs > 1)
                    geno_sub[:, af_is_not_maf] = 2 - geno_sub[:, af_is_not_maf]

                    E = env_data[:, j]
                    geno_EDGE2 = transform_genotypes(geno_sub, p, E)
                    E_perm = np.random.choice(E, len(E), replace = False)
                    geno_EDGE2_null = transform_genotypes(geno_sub, p, E_perm)
                    #geno_EDGE22 = transform_genotypes(geno_sub, p, E, False)
                    for g1, g2, g2_null, rsID in zip(geno_sub.T, geno_EDGE2.T,  geno_EDGE2_null.T, rsIDs_j):
                        EDGE2_geno_pheno_sets[i].append(g2)
                        geno_pheno_sets[i].append(g1)
                        val_indices = np.logical_or(np.isnan(g1), np.isnan(p2)) == False
                        X = np.array([g1, E, g1*E, np.ones(len(E))]).T
                        X_sub = np.array([g1, E, np.ones(len(E))]).T
                        model = sm.OLS(p[val_indices], X[val_indices]).fit()
                        model_sub = sm.OLS(p[val_indices], X_sub[val_indices]).fit()
                        linear_p_val = model.compare_lr_test(model_sub)[1]
                        standard_p_values_GxE.append(linear_p_val)
                        standard_p_values_E_ind.append(all_E_names[j])
                        standard_p_values_rsids.append(rsID)
                        
                        if rsID == "rs13118211":
                            X2 = np.array([g2, E, np.ones(len(E))]).T
                            model2 = sm.OLS(p, X2).fit()

                            m1, b1 = linregress(E, p)[0:2]
                            m2, b2 = linregress(E, g2)[0:2]
                            p_res = p - (E*m1 + b1)
                            g2_res = g2 - (E*m2 + b2)
                            r, void = pearsonr(g2_res, p_res)

                            m1_null, b1_null = linregress(E_perm, p)[0:2]
                            m2_null, b2_null = linregress(E_perm, g2_null)[0:2]
                            p_res_null = p - (E_perm*m1_null + b1_null)
                            g2_res_null = g2_null - (E_perm*m2_null + b2_null)
                            r_null, void = pearsonr(g2_res_null, p_res_null)

                            N = int(len(g2_res)/3)
                            sorted_indices = np.argsort(g2_res)
                            X_vals = np.unique(g2_res)
                            X_sorted = g2_res[sorted_indices]
                            X_edge = [np.mean(X_sorted[i*N:(i+1)*N]) for i in range(2)] + [np.mean(X_sorted[2*N:])]

                            Y_sorted = p_res[sorted_indices]
                            Y_edge_dists = [Y_sorted[i*N:(i+1)*N] for i in range(2)] + [Y_sorted[2*N:]]
                            Y_edge_mean = [np.mean(dist) for dist in Y_edge_dists]
                            Y_edge_median = [np.median(dist) for dist in Y_edge_dists]

                            Y_edge_special = [np.mean(np.percentile(dist, [10, 15, 20, 80, 85, 90])) for dist in Y_edge_dists]
                            plt.plot(X_edge, Y_edge_mean, "*", label = "means")
                            plt.plot(X_edge, Y_edge_median, "*", label = "medians")
                            plt.plot(X_edge, Y_edge_special, "*", label = "average of percentiles 10, 15, 20, 80, 85, 90")
                            plt.legend()
                            plt.savefig("step10.png")
                            plt.clf()

                            N = int(len(g2_res_null)/3)
                            sorted_indices_null = np.argsort(g2_res_null)
                            X_vals_null = np.unique(g2_res_null)
                            X_sorted_null = g2_res_null[sorted_indices_null]
                            X_edge_null = [np.mean(X_sorted_null[i*N:(i+1)*N]) for i in range(2)] + [np.mean(X_sorted_null[2*N:])]

                            Y_sorted_null = p_res_null[sorted_indices_null]
                            Y_edge_dists_null = [Y_sorted_null[i*N:(i+1)*N] for i in range(2)] + [Y_sorted_null[2*N:]]
                            Y_edge_mean_null = [np.mean(dist) for dist in Y_edge_dists_null]
                            Y_edge_median_null = [np.median(dist) for dist in Y_edge_dists_null]

                            figure, axes = plt.subplots(nrows = 2)
                            axes[0].plot(X_edge, Y_edge_mean, "*", label = "EDGE2 effect (r = " + str(r)[0:6] + ")")
                            axes[0].plot(X_edge_null, Y_edge_mean_null, "*", label = "permuted EDGE2 effect (r = " + str(r_null)[0:6] + ")")
                            axes[0].set_xlabel("mean encoding for bottom/middle/top\n thirds of the EDGE2 encodings")
                            axes[0].set_ylabel("mean phenotype\n per genotype")
                            axes[0].legend()  
                            axes[0].title.set_text("EDGE2 encoding (p = 6.62E-14)")

                            E_halves = [E <= np.median(E), E > np.median(E)]
                            labels = ["no smoking", "smoking"]
                            for E_ind, E_label in zip(E_halves, labels):
                                ge = g1[E_ind]
                                Ye = p[E_ind]
                                X = [0,1,2]
                                Ye_means = [np.mean(Ye[ge == i]) for i in X] 
                                axes[1].plot(X, Ye_means, "*", label = E_label)
                            axes[1].set_xlabel("additively encoded genotype")
                            axes[1].set_ylabel("mean phenotype per\n genotype and environment")
                            axes[1].legend()
                            axes[1].title.set_text("additive encoding (p = " + str(linear_p_val)[0:4] + ")")
                            plt.xticks([0,1,2], [0,1,2])
                            figure.tight_layout(pad = 2.1)
                            figure.set_figheight(6.5)
                            figure.set_figwidth(5)
                            plt.savefig("rough_figure.png")
                            plt.clf()

                            for i, dist in enumerate(Y_edge_dists):
                                plt.hist(dist, bins = 100, label = "distribution " + str(i + 1))
                                plt.savefig("step10_dist" + str(i + 1) + ".png")
                                plt.clf()
                               

temp_df = [standard_p_values_rsids, standard_p_values_E_ind, standard_p_values_GxE]
temp_df = pd.DataFrame(np.array(temp_df).T)
temp_df.columns = ["rsID", "env_factor", "pval_null"]
GxE_info = GxE_effects.merge(temp_df, on = ["rsID", "env_factor"], how = "inner")
p_alt = -np.log10(GxE_info["pEDGE2"].to_numpy())
p_null = -np.log10(GxE_info["pval_null"].to_numpy(dtype = float))
env_factors = GxE_info["env_factor"].to_numpy()
sig_ind = (p_null >= -np.log10(pb))

Interesting_SNPs = GxE_info[sig_ind == False]
Interesting_SNPs.index = np.arange(len(Interesting_SNPs))
diffs = -np.log10(Interesting_SNPs["pEDGE2"].to_numpy(dtype = float)) + np.log10(Interesting_SNPs["pval_null"].to_numpy(dtype = "float"))
best_SNP_ind = np.argsort(diffs)
Interesting_SNPs.loc[best_SNP_ind[-2:], :]

'''
# some commentary
# p_null is the GxE p value for a linear interaction test
# p_alt is the permutation test that I use. 
# sig_ind are the indices where p_null is significant
p_alt1, p_alt2 = p_alt[sig_ind], p_alt[sig_ind == False]
p_null1, p_null2 = p_null[sig_ind], p_null[sig_ind == False]
p_total = wilcoxon(p_alt - p_null, alternative = 'greater')

# As expected, p_alt greatly outperforms p_null when p_null is not significant
p_not_sig_only = wilcoxon(p_alt2 - p_null2, alternative = 'greater')

# As expected, p_alt does not outperform p_null when p_null is significant
p_sig_only = wilcoxon(p_alt1 - p_null1, alternative = 'greater')

# We expect p_alt to outperform p_null more when pnull is insignificant, and this is the case
method_diff = mwu(p_alt2 - p_null2, p_alt1 - p_null1, alternative = 'greater')

# We expect p_alt1 and p_alt2 not to be significantly different, and they are not.
# This is because p_alt's significance should be uneffected by p_null's significance. 
p_alt_diff = mwu(p_alt1, p_alt2, alternative = 'two-sided')

spearmanr(p_alt1, p_null1)
spearmanr(p_alt2, p_null2)
spearmanr(p_alt, p_null)
'''

r2_vals = []
p_vals = []
PRS_scores = []
r2_vals_null = []
p_vals_null = []
#for i in tqdm(phenotypes[:-1]):
for i in tqdm(phenotypes):
    N = len(env_data)
    ph = COPY(QTL_phenotypes[:, i])
    ph[np.isnan(ph)] = np.nanmedian(ph)
    if len(EDGE2_geno_pheno_sets[i]) > 0:
        geno_EDGE2 = np.concatenate([np.array(EDGE2_geno_pheno_sets[i]).T, env_data, np.ones((N, 1))], axis = 1)
        geno_null =  np.concatenate([np.array(geno_pheno_sets[i]).T, np.ones((N, 1))], axis = 1)
        for k in range(len(geno_null[0])): geno_null[np.isnan(geno_null[:, k]), k] = np.nanmedian(geno_null[:, k])
        model1 = sm.OLS(ph, geno_EDGE2).fit()
        L = len(env_data[0]) + 1
        beta = model1.params[:-L]
        ph_PRS = np.matmul(geno_EDGE2[:, :-L], beta)
        r, p = pearsonr(ph, ph_PRS)
        PRS_scores.append(ph_PRS)
        r2_vals.append(r**2)
        p_vals.append(p)
        model2 = sm.OLS(ph, geno_null).fit()
        r2_vals_null.append(model2.rsquared) 

path = "../step7_adjust_HF_for_covariates_NN/y.txt"
data = pd.read_csv(path, delimiter = "\t") 
data = data.merge(eids, on = "eid", how = "inner")
y = data['any_HF'].to_numpy(dtype = int)
X = np.array(PRS_scores).T
#X = []
#for i in EDGE2_geno_pheno_sets: X += i
#X = np.array(X).T
case_inds = np.where(y == 1)[0]
control_inds = np.where(y == 0)[0]
keepers1 = np.repeat(case_inds, 2)
keepers2 = np.random.choice(control_inds, size = len(case_inds), replace = False)
keepers = np.concatenate([keepers1, keepers2])
pd.DataFrame(X).to_csv("step10e_NN_PRS.txt", sep = "\t", header = False, index = False)
X2, y2 = X[keepers], y[keepers]

n_est = [10, 20, 30, 40, 50]
lr = [0.1, 0.2, 0.3, 0.4, 0.5] 
max_depth = [1, 2, 3, 4, 5]
combinations = []
mean_r2_vals = []
for i in tqdm(n_est):
    for j in lr:
        for k in max_depth:
            combinations.append((i, j, k))
            clf = GBC(n_estimators = i, learning_rate = j, max_depth = k)
            r2 = cross_validate(X2, y2, 10, clf)
            mean_r2_vals.append(np.mean(r2))

step10e_params = pd.DataFrame(np.array(combinations[np.argmax(mean_r2_vals)]).T)
step10e_params.to_csv("step10e_params.txt", sep = "\t", header = True, index = False)

p1, p2, p3 = combinations[np.argmax(mean_r2_vals)]
clf1 = GBC(n_estimators = p1, learning_rate = p2, max_depth = p3)
clf2 = LogisticRegression(C = 1E20, solver = 'sag')
r21, r22 = cross_validate(X2, y2, 30, clf1, clf2)
pval = wilcoxon(r21 - r22, alternative = 'greater')[1]

step10e_r_vals = pd.DataFrame(np.array([r21, r22]).T)
step10e_r_vals.columns = ["GBC", "LR"]
step10e_r_vals.to_csv("step10e_r_vals.txt", sep = "\t", header = True, index = False)

pdb.set_trace()

prefix = "../step10_get_significant_SNPs_PCA/"
X2_PCA = pd.read_csv(prefix + "step10e_PCA_PRS.txt", delimiter = "\t", header = None).to_numpy()
X3 = np.concatenate([X2, X2_PCA[keepers, :]], axis = 1)
r31, r32 = cross_validate(X3, y2, 30, clf1, clf2)
pval = wilcoxon(r31 - r32, alternative = 'greater')[1]

step10e_r_vals = pd.DataFrame(np.array([r31, r32]).T)
step10e_r_vals.columns = ["GBC", "LR"]
step10e_r_vals.to_csv("step10e_r_vals_expanded.txt", sep = "\t", header = True, index = False)

pdb.set_trace()

family = sm.genmod.families.Binomial()
family.link = sm.genmod.families.links.logit()
X4 = np.concatenate([X2, np.ones((len(X2), 1))], axis = 1)
model = sm.GLM(y2, X4, family = family)
model_results = model.fit()
sig_indices = np.arange(len(X4[0]))[model_results.pvalues < 0.05]
model2 = sm.GLM(y2, X4[:, sig_indices], family = family)
model2_results = model2.fit()
LR_test_statistic = 2*(model_results.llf - model2_results.llf)
model_diff_p = chi2.sf(LR_test_statistic, len(X4[0]) - len(sig_indices))
if np.any(model2_results.pvalues > 0.05) or model_diff_p < 0.05:
    pdb.set_trace()

component_inds = [i + 1 for i in sig_indices]
component_inds[-1] = "intercept"
step10e_sig_components = pd.DataFrame(np.array([component_inds, model2_results.pvalues, model2_results.params], dtype = object).T)
step10e_sig_components.columns = ["component", "p value", "beta"]
step10e_sig_components = step10e_sig_components.sort_values("p value")
step10e_sig_components.to_csv("step10e_sig_components.txt", sep = "\t", header = True, index = False)

family = sm.genmod.families.Binomial()
family.link = sm.genmod.families.links.logit()
X5 = np.concatenate([X3, np.ones((len(X3), 1))], axis = 1)
model = sm.GLM(y2, X5, family = family)
model_results = model.fit()
sig_indices = np.arange(len(X5[0]))[model_results.pvalues < 0.05]
model2 = sm.GLM(y2, X5[:, sig_indices], family = family)
model2_results = model2.fit()
LR_test_statistic = 2*(model_results.llf - model2_results.llf)
model_diff_p = chi2.sf(LR_test_statistic, len(X5[0]) - len(sig_indices))
if np.any(model3_results.pvalues > 0.05) or model_diff_p < 0.05:
    pdb.set_trace()

# TODO: put PCA vs logistic PCA labels back
component_inds = [i + 1 for i in sig_inds]
component_inds[-1] = "intercept"
step10e_sig_components = pd.DataFrame(np.array([component_inds, model2_results.pvalues, model2_results.params], dtype = object).T)
step10e_sig_components.columns = ["component", "p value", "beta"]
step10e_sig_components = step10e_sig_components.sort_values("p value")
step10e_sig_components.to_csv("step10e_sig_components_expanded.txt", sep = "\t", header = True, index = False)

'''
colnames = ["mean EDGE2 -logp value nonlinear effect", "mean standard -logp value nonlinear effect", "diff p value"]
colnames += ["mean EDGE2 -logp value linear effect", "mean standard -logp value linear effect", "diff p value"]
colnames += ["mean EDGE2 -logp value nonlinear effect", "mean EDGE2 -logp value linear effect", "diff p value"]
colvals = [np.mean(p_alt2), np.mean(p_null2), p_not_sig_only[1]]
colvals += [np.mean(p_alt1), np.mean(p_null1), p_sig_only[1]]
colvals += [np.mean(p_alt1), np.mean(p_alt2), p_alt_diff[1]]
df = pd.DataFrame(np.array([colnames, colvals], dtype = object).T)
df.to_csv("step10e_global_results.txt", sep = "\t", header = False, index = False)
'''

GxE_info["logp diff"] = np.log10(GxE_info["pval"].astype(float)) - np.log10(GxE_info["pval_null"].astype(float))
GxE_info = GxE_info.sort_values("logp diff")
GxE_info.to_csv("step10e_SNP_info.txt", sep = "\t", header = False, index = False)
