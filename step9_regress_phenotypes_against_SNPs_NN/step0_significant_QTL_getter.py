import numpy as np
import pandas as pd
import argparse
import os
from functools import reduce
from copy import deepcopy as COPY
from bed_reader import open_bed
import statsmodels.api as sm
from scipy.stats import pearsonr
from scipy.stats import linregress
from scipy.stats import mode
from scipy.stats import t
from scipy.stats import chi2
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
import pdb
from tqdm import tqdm

def nanlinregress(x, y):
    xval_indices = (np.isnan(x) == False)
    yval_indices = (np.isnan(y) == False)
    val_indices = np.logical_and(xval_indices, yval_indices)
    x, y = x[val_indices], y[val_indices]
    return(linregress(x, y))

# taken from https://online.stat.psu.edu/stat462/node/131/
# and https://online.stat.psu.edu/stat462/node/137/
def np_OLS(X_full, X_less, y):
  
    Betas2_full, SSE_full = np.linalg.lstsq(X_full, y, rcond = None)[0:2]
    Betas2_less, SSE_less = np.linalg.lstsq(X_less, y, rcond = None)[0:2]
    d_full = len(y) - 3
    T = ((SSE_less - SSE_full)/(SSE_full/d_full))**(1/2)
    p_val1 = 2*t(d_full).sf(T)[0]

    '''
    python step0_significant_QTL_getter.py --chr X --pheno 0
    # just confirms my p values are computed correctly
    model = sm.OLS(y, X_full)
    p_val2 = model.fit().pvalues[0]
    if not np.round(p_val1, 6) == np.round(p_val2, 6):
        pdb.set_trace()
    else:
        print(p_val1)
        print(p_val2)
        print("\n\n")
    '''
    return(Betas2_full, p_val1)


if not os.path.exists("QTL_output"):
    try: 
        os.mkdir("QTL_output")
    except:
        pass

parser = argparse.ArgumentParser()
parser.add_argument('--chr', nargs = 1, type = str, action = "store", dest = "chr")
parser.add_argument('--pheno', nargs = 1, type = int, action = "store", dest = "pheno")
args = parser.parse_args()
chr = (args.chr)[0]
pheno = (args.pheno)[0]

filename = "QTL_output/QTL_effects_chr" + chr + "_P" + str(pheno) + ".txt"
if os.path.exists(filename):
    exit()

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
new_fam_path_prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr"
new_fam_paths = [new_fam_path_prefix + i + ".fam" for i in chromosomes]
new_fams = [pd.read_csv(path, delim_whitespace = True, header = None)[1].to_numpy() for path in new_fam_paths]
main_path = np.array(new_fam_paths)[np.array(chromosomes) == chr][0]
new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
new_fam_intersect = reduce(np.intersect1d, new_fams)
is_in_intersect = np.isin(new_fam_main, new_fam_intersect)
sorted_main_indices = np.argsort(new_fam_main)
sorted_indices = sorted_main_indices[is_in_intersect[sorted_main_indices]]
new_fam = new_fam_main[sorted_indices]

#                                                                   ex             ex               ex          gender
# maybe start with these: ['pack-years', 'annual-consumption', '874-average',  '894-average', '914-average', '22001-0.0']
# env_name = 'pack-years'
# env_name = 'annual-consumption'
# env_name = ['874-average', '894-average', '914-average']
env_name = '22001-0.0'
path = "../step7_adjust_HF_for_covariates_NN/env_factors_for_step9.txt"
gender_path = "../step7_adjust_HF_for_covariates_NN/X.txt"
gender_data = pd.read_csv(gender_path, delimiter = "\t", usecols = ["eid", "22001-0.0"], dtype = int)
if env_name == ['874-average', '894-average', '914-average']:
    env_data = pd.read_csv(path, delimiter = "\t", usecols = (["eid"] + env_name))
    exercise = env_data[env_name].to_numpy()
    pca = PCA(n_components = 1)
    pca.fit(exercise)
    exercise_score = pca.transform(exercise).reshape(-1)
    env_data["exercise_score"] = exercise_score
    env_data = env_data[["eid", "exercise_score"]]
    env_name = "exercise_score"
else:
    env_data = pd.read_csv(path, delimiter = "\t", usecols = ["eid", env_name])
env_data["eid"] = env_data.loc[:, "eid"].astype(int)
new_fam_df = pd.DataFrame(np.array([new_fam, np.arange(len(new_fam))]).T)
new_fam_df.columns = ["eid", "index"] 
gender_data = gender_data.merge(new_fam_df, on = "eid", how = "inner").sort_values(by = "index")
env_data = env_data.merge(new_fam_df, on = "eid", how = "inner").sort_values(by = "index")
is_male = gender_data["22001-0.0"].to_numpy(dtype = "bool")
env_factor = env_data[env_name].to_numpy()

old_fam_path = "../step4_remove_relatives/UKB_samples_unrelated.fam"
old_fam = pd.read_csv(old_fam_path, delim_whitespace = True, header = None)[0].to_numpy()
old_in_new = np.isin(old_fam, new_fam)
# changing the definition of heart failure led to a slight change in the eids, which cannot be reflected in step 8
# only unrelated individuals under the previous procedure (step 8) and the current procedure (step 4) are included
# This does include all all-cause heart failure cases and only removes a few thousand controls.
new_in_old = np.isin(new_fam, old_fam)
if not np.all(old_fam == np.sort(old_fam)):
    print("error1: code expects the fam file from step 4 to have sorted eids")
    exit()
if not np.all(new_fam[new_in_old] == old_fam[old_in_new]):
    print("error2: code expects the fam file from step 4 to have sorted eids")
    exit()

geno_path = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr" + chr
pheno_path = "../step7_adjust_HF_for_covariates_NN/phenotypes_for_step9.txt"
rsIDs = pd.read_csv(geno_path  + ".bim", delim_whitespace = True, header = None)[1]
phenotypes_info = pd.read_csv(pheno_path, delimiter = "\t").loc[old_in_new, :]
phenotypes_info[["eid"]].to_csv("QTL_phenotypes_eids.txt", sep = "\t", header = False, index = False)
phenotypes = phenotypes_info.loc[:, np.arange(16).astype(str)].to_numpy()
if chr == '1' and pheno == 0:
    pd.DataFrame(phenotypes).to_csv("QTL_phenotypes.txt", sep = "\t", header = False, index = False)
    pd.DataFrame(is_male).to_csv("is_male.txt", sep = "\t", header = False, index = False)
    # env_factors = pd.read_csv(path, delimiter = "\t")
    # env_factors = env_factors.loc[old_in_new, :]
    # env_factors.to_csv("env_factors_cleaned.txt", sep = "\t", header = True, index = False)
phenotypes = phenotypes[:, pheno]
if chr == "Y":
    phenotypes = phenotypes[is_male]

prefix = "genotype_metadata/genotype_metadata_chr"
cols = ["SNP", "C(HOM A1)", "C(HET)", "C(HOM A2)", "C(MISSING)"]
all_missingness = pd.read_csv(prefix + chr + ".frqx", delimiter = "\t", usecols = ["SNP", "C(MISSING)"])
all_afs = pd.read_csv(prefix + chr + ".frq", delim_whitespace = True, usecols = ["SNP", "MAF"])
if not np.all(all_missingness["SNP"] == rsIDs) and np.all(all_afs["SNP"] == rsIDs):
    print("exiting: frq and frqx rsID are not ordered as expected")
all_missingness = all_missingness["C(MISSING)"].to_numpy()/len(new_fam_main)
all_afs = all_afs["MAF"].to_numpy()
all_mafs = maf2 = np.min(np.array([all_afs, 1 - all_afs]), axis = 0)
non_missing_indices = np.arange(len(rsIDs))[all_missingness <= 0.02]
non_constant_indices =  np.arange(len(rsIDs))[all_mafs >= 0.005]
good_indices = np.intersect1d(non_missing_indices, non_constant_indices)
intervals = np.cumsum([0] + int(len(good_indices)/100)*[100] + [len(good_indices)%100])


N = len(good_indices)
new_rsIDs = []
out = np.zeros((N, 8))
i_final = 0
for j in tqdm(range(len(intervals) - 1)):

    col_indices = good_indices[intervals[j]:intervals[j + 1]]
    rsIDs_chunk = rsIDs[col_indices].to_numpy()
    new_rsIDs += rsIDs_chunk.tolist()

    # my_rsID = 'rs12687942'
    # if my_rsID in rsIDs_chunk:
    #     pdb.set_trace()
    #     genotypes = open_bed(geno_path  + ".bed", num_threads = 1).read(np.s_[sorted_indices, col_indices])
    #     genotypes = genotypes[new_in_old, :]
    #     ind = np.where(rsIDs_chunk == my_rsID)[0][0]
    #     zz = nanlinregress(genotypes[:, ind], phenotypes)
    # else:
    #     continue

    genotypes = open_bed(geno_path  + ".bed", num_threads = 1).read(np.s_[sorted_indices, col_indices])
    # removes rows from step 8 not in modified step4
    genotypes = genotypes[new_in_old, :]
    is_male2 = is_male
    env_factor2 = env_factor
    if chr == "Y":
        genotypes, is_male2, env_factor2 = genotypes[is_male], is_male[is_male], env_factor[is_male]

    afs = all_afs[col_indices]
    mafs = all_mafs[col_indices]
    af_is_not_maf = (afs != mafs)
    genotypes[:, af_is_not_maf] = 2 - genotypes[:, af_is_not_maf]
    valued_indices = np.logical_or(np.isnan(genotypes), np.isnan(phenotypes).reshape(-1,1)) == False

    env_is_not_gender = np.any(env_factor2 != is_male2)
    if chr == "X" and env_is_not_gender:
        all_genotypes = [genotypes[is_male2], genotypes[is_male2 == False]]
        all_phenotypes = [phenotypes[is_male2], phenotypes[is_male2 == False]]
        all_env = [env_factor2[is_male2], env_factor2[is_male2 == False]]
        all_valued_indices = [valued_indices[is_male2], valued_indices[is_male2 == False]]
        all_vals = zip(all_genotypes, all_phenotypes, all_env, all_valued_indices, [0, 1])
        L = len(genotypes[0])
        p_main, p_null1, p_null2 = np.zeros((2, L)), np.zeros((2, L)), np.zeros((2, L))
        num_groups = 2
    else:
        all_vals = zip([genotypes], [phenotypes], [env_factor2], [valued_indices], [0])
        num_groups = 1
    for geno, pheno_vals, env, val_inds, ind in all_vals:
        i = COPY(i_final)
        for k in range(len(geno[0])):
            g = geno[val_inds[:, k], k]
            out[i, 7] = len(g)
            p = pheno_vals[val_inds[:, k]]
            E1 = env[val_inds[:, k]]
            g2_vals = (g == 2)
            is_standard = np.sum(g2_vals, axis = 0) > 1000 and chr != "Y"
            if is_standard: 

                # if rsIDs_chunk[k] == my_rsID:
                #      pdb.set_trace()

                g1_vals = (g == 1)
                X = np.array([g1_vals, g2_vals, (g1_vals)*E1, (g2_vals)*E1, E1, np.ones(len(g1_vals))]).T
                X2 =[g1_vals, g2_vals]
                if chr in ["X", "XY"] and env_is_not_gender == False:
                    X2 += [is_male2[val_inds[:, k]]]
                X2 = np.array(X2 + [np.ones(len(g1_vals))]).T
                y = (p - np.mean(p))

                Betas = np.linalg.lstsq(X, y, rcond = None)[0][:-1]
                enc1, enc2 = Betas[0] + Betas[2]*E1[g1_vals], Betas[1] + Betas[3]*E1[g2_vals] 
                g_old = COPY(g)
                g[g1_vals] = enc1     
                g[g2_vals] = enc2   

                Betas2 = np.linalg.lstsq(X2, y, rcond = None)[0][:-1]
                g2 = COPY(g_old)
                g2[g1_vals] = Betas2[0]  
                g2[g2_vals] = Betas2[1]

                ge = [g, E1]
                ge = np.array(ge + [np.ones(len(g))]).T
                params, p_val = np_OLS(ge, ge[:, 1:], p)
                if num_groups == 1:
                    out[i, 1] = p_val
                else:
                    p_main[ind][k] = p_val

                # if rsIDs_chunk[k] == my_rsID:
                #     pdb.set_trace()

            elif not (chr == "Y" and len(np.unique(E1)) == 1):
                g[g2_vals] = 1
                g1_vals = (g == 1)
                X = np.array([g1_vals, (g1_vals)*E1, E1, np.ones(len(g1_vals))]).T
                y = (p - np.mean(p))

                Betas = np.linalg.lstsq(X, y, rcond = None)[0][:-1]
                enc = Betas[0] + Betas[1]*E1[g1_vals]   
                g_old = COPY(g)
                g[g1_vals] = enc 

                ge = np.array([g, E1, np.ones(len(g))]).T
                params, p_val = np_OLS(ge, ge[:, 1:], p)
                if num_groups == 1:
                    out[i, 1] = p_val
                else:
                    p_main[ind][k] = p_val


            else:
                g_old = COPY(g)
                if num_groups == 1:
                    out[i, 1] = linregress(g, p)[3]
                else:
                    p_main[ind][k] = linregress(g, p)[3]


            if is_standard:
                if num_groups == 1:
                    out[i, 2] = linregress(g2, p)[3]
                else:
                    p_null1[ind][k] = linregress(g2, p)[3]
            else:
                out[i, 2] = -1
            g_norm = (g_old - np.mean(g_old))/np.std(g_old)
            slope, void, void, null2_p_val = linregress(g_norm, p)[0:4]

            out[i, 3 + ind] = null2_p_val
            out[i, 5 + ind] = slope
            i += 1
    if num_groups == 2:
        out[(i - k - 1):(i), 1] = chi2(4).sf(-2*np.sum(np.log(p_main), axis = 0))
        if is_standard:
            out[(i - k - 1):(i), 2] = chi2(4).sf(-2*np.sum(np.log(p_null1), axis = 0))
    i_final = COPY(i)

#pdb.set_trace()
o_df = pd.DataFrame(out)
# if chromosome is X, then null2 is male and null2_alt is female. Otherwise, it's all samples. 
o_df.columns = ["rsID", "p_main", "p_null1", "p_null2", "p_null2_alt", "p_null2_beta", "p_null2_alt_beta", "N"]
o_df["rsID"] = new_rsIDs
o_df.to_csv("QTL_output/QTL_effects_chr" + chr + "_P" + str(pheno) + ".txt", sep = "\t", header = True, index = False)