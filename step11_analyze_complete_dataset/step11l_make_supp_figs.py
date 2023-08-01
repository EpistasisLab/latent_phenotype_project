import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import rankdata
from copy import deepcopy as COPY
import os 
import pdb

# note that step11g_make_heritability_figure.py makes table S4.

# Table S1 starts here

fname = "../step7_adjust_HF_for_covariates_PCA/y.txt"
Y_df = pd.read_csv(fname, delimiter = "\t", header = 0)
all_cols = Y_df.columns.to_numpy()
ICD_cols = all_cols[np.isin(all_cols, ["eid", 'any_HF']) == False]
ICD_codes = np.array([code.split("_")[1] for code in ICD_cols])
ICD_codes = np.array([code[:3] + "." + code[3:] if len(code) > 3 else code for code in ICD_codes])
ICD_categories = np.array([code[0:2] if code[0] == "I" else code[0] for code in ICD_codes])
ICD_descriptions = ["rheumatic illness", "hypertensive diseases", "ischaemic/embolic heart diseases"]
ICD_descriptions += ["pericardial and valve diseases", "conduction disorders, cardiomyopathies, and arrhythmias", "heart failure"]
ICD_descriptions += ["haemorrhages, stenosis, occlusions, infarctions", "atherosclerosis, aneurysms, embolisms", "thrombosis, varicose veins, haemorrhoids"]
ICD_descriptions += ["hypotension", "diseases of the respiratory system", "miscellaneous abnormal clinical findings (correlated to AHF)"]
codes, categories, subcategories = [], [], []
for i in np.unique(ICD_categories):
    inds = (ICD_categories == i)
    codes.append(", ".join(ICD_codes[inds]))
    categories.append(i[0])
    if len(i) > 1: subcategories.append(i[1])
    else: subcategories.append("all")

df = pd.DataFrame(np.array([categories, subcategories, ICD_descriptions, codes]).T)
df.columns = ["ICD10 code category", "sub-category", "description", "ICD10 codes included"]
df.to_csv("tableS1.txt", sep = "\t", header = True, index = False)

# Figure S1 starts here
prefixes = ["/home/greggj/pleiotropy_and_GxE/step9_regress_phenotypes_against_SNPs_PCA/QQ_plots_smoking/",
           "/home/greggj/pleiotropy_and_GxE/step9_regress_phenotypes_against_SNPs_logistic_PCA/QQ_plots_smoking/",
           "/home/greggj/pleiotropy_and_GxE/step9_regress_phenotypes_against_SNPs_NN/QQ_plots_smoking/"]

IF_filenames = []
for prefix in prefixes: IF_filenames += [prefix + name for name in os.listdir(prefix) if "IF" in name]
IFs = np.array([pd.read_csv(name, delimiter = "\t", header = None).loc[0, 0] for name in IF_filenames])
IFs = np.array([pd.read_csv(name, delimiter = "\t", header = None).loc[0, 0] for name in IF_filenames])
plt.figure(figsize=(24, 16))
plt.hist(IFs, bins = 48)
plt.title("Inflation Factors for\nLatent Phenotypes\' Main Effects", fontsize = 64)
plt.xlabel('Inflation Factor', fontsize=40)
plt.ylabel('Bin Count', fontsize=40)
plt.xlim([0.85, 1.55])
plt.xticks(fontsize=32)
plt.yticks(fontsize=32)
plt.subplots_adjust(top=0.84)
plt.savefig("figureS1.png")
plt.clf()

# Table S2 starts here
prefix = "/home/greggj/pleiotropy_and_GxE/step7_adjust_HF_for_covariates_PCA"
imputation_helpers = pd.read_csv(prefix + "/step7d_imputation_helper_cols.txt", delimiter = "\t")
is_imputed_prefix = ""
is_imputed_suffix = ""
features = []
descriptions = []
for i in imputation_helpers["0"].to_numpy():
    possible_field = ''.join(i.split("-0.0"))
    if "is_imp" in i:
        is_imputed_suffix = " is imputed"
        is_imputed_prefix = "whether or not a "
        possible_field = possible_field.split("_is_")[0]
    try:
        field = int(possible_field)
        features.append(possible_field + is_imputed_suffix)
        descriptions.append(is_imputed_prefix + "value of UK biobank field " + possible_field + is_imputed_suffix)
    except:
        features.append(possible_field + is_imputed_suffix)
        if len(possible_field) <= 5:
            descriptions.append("whether or not patient has UK biobank ICD10 code " + possible_field)
        elif len(possible_field.split("-")) == 2:
            part1, part2 = possible_field.split("-")
            descriptions.append(is_imputed_prefix + "value of UK biobank field " + part1 + " does or does not equal " + part2.split(".")[0] + is_imputed_suffix)
        else:
            descriptions.append(is_imputed_prefix + "value of custom feature " + possible_field + is_imputed_suffix)
    is_imputed_prefix = ""
    is_imputed_suffix = ""

tableS2 = pd.DataFrame(np.array([features, descriptions]).T)
tableS2.columns = ["feature name", "feature description"]
tableS2.to_csv("tableS2.txt", sep = "\t", header = True, index = False)

# Table S3 starts here

prefixes = ["/home/greggj/pleiotropy_and_GxE/step7_adjust_HF_for_covariates_PCA/"]
prefixes += ["/home/greggj/pleiotropy_and_GxE/step7_adjust_HF_for_covariates_logistic_PCA/"]
prefixes += ["/home/greggj/pleiotropy_and_GxE/step7_adjust_HF_for_covariates_NN/"]

colnames = ["simulated missingness pattern", "pack-years of smoking", "annual alcohol consumption"]
colnames += [ "UKB field 874", "UKB field 894", "UKB field 914"]

patterns = ["random values", "highest values", "lowest values"]
patterns += ["highest first principle component values", "lowest first principle component values"]
for prefix, letter in zip(prefixes, ["a", "b", "c"]):
    
    table = pd.read_csv(prefix + "corr_sets_0.05neighbors.txt", delimiter = "\t")
    table.columns = colnames
    table["simulated missingness pattern"] = patterns
    table.to_csv("tableS3" + letter + ".txt", sep = "\t", header = True, index = False)

# Figure S2 starts here

base = "../step10_get_significant_SNPs_"
all_files  = []
for model, letter in zip(["PCA", "logistic_PCA", "NN"], ["a", "b", "c"]):
    prefix = base + model + "/"
    QTL_suffixes = ["alcohol", "exercise", "smoking", "gender"]
    paths = [prefix + "hits_QTL_" + i for i in QTL_suffixes]
    files = []
    for path in paths:
        env_factor = path.split("_")[-1]
        all_paths = [path + "/" + i for i in os.listdir(path)]
        file = pd.concat([pd.read_csv(path, delimiter = "\t") for path in all_paths])
        file["env_factor"] = env_factor
        files.append(file)
    QTL_files = pd.concat(files)[["rsID", "pheno_index", "p_main", "p_null2", "env_factor"]]
    QTL_files["p"] = (QTL_files["p_main"]/QTL_files["p_null2"]).to_numpy()
    QTL_files = QTL_files[["rsID", "pheno_index", "p", "env_factor"]]

    GxE_effects_df = pd.read_csv(prefix + "step10e_GxE_effects_df.txt", delimiter = "\t")
    all_files = QTL_files.merge(GxE_effects_df, on = ["rsID", "pheno_index", "env_factor"], how = "inner")
    plt.figure(figsize=(10, 7.5))
    pvals = all_files[["p", "pEDGE2"]]
    corr, void = pearsonr(-np.log10(pvals["p"]), -np.log10(pvals["pEDGE2"]))
    outlier = np.argmin(pvals["p"])
    inds = np.setdiff1d(np.arange(len(pvals)), [outlier])
    plt.plot(-np.log10(pvals["p"])[inds], -np.log10(pvals["pEDGE2"])[inds], "ok", label = "r = " + str(np.round(corr, 2)))
    plt.legend(fontsize=20, loc = "upper left")
    titlemodel = COPY(model)
    if model == "logistic_PCA": titlemodel = "logistic PCA"
    plt.title("permutation GxE p-value vs proxy (" + titlemodel + ")", fontsize = 24, y=1.045)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.subplots_adjust(top=0.88)
    plt.xlabel("-log10(proxy p-value)", fontsize=20)
    plt.ylabel("-log10(permutation p-value)", fontsize=20)
    plt.savefig("figureS2" + letter + ".png")
    plt.clf()

# Table S4 starts here

maf_base = "../step9_regress_phenotypes_against_SNPs_PCA/genotype_metadata"
maf_fnames = [maf_base + "/" + i for i in os.listdir(maf_base) if ".frq" in i and ".frqx" not in i]
mafs = pd.concat([pd.read_csv(name, delim_whitespace = True, usecols = ["SNP", "MAF"]) for name in maf_fnames])
mafs["MAF"] = np.min(np.array([mafs["MAF"], 1 - mafs["MAF"]]), axis = 0)
mafs.columns = ["rsID", "minor allele frequency"]

base = "../step10_get_significant_SNPs_"
n, letters = 0, ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o"] 
for model in ["PCA", "logistic_PCA", "NN"]:
    for E in ["main", "GxSmoking", "GxGender", "GxAlcohol", "GxExercise"]:
        
        if E == "main": E2 = "GxSmoking"
        else: E2 = E
        pos_base = "../step10_get_significant_SNPs_" + model + "/significant_SNPs_plink_files_" + E2.split("Gx")[-1].lower()
        pos_fnames = [pos_base + "/" + i for i in os.listdir(pos_base) if ".bim" in i]
        pos = pd.concat([pd.read_csv(name, delim_whitespace = True, header = None) for name in pos_fnames])[[1, 0, 3]]
        pos.columns = ["rsID", "chromosome", "position"]

        fname = base + model + "/" + "rsIDs_" + E + "_effects.txt"
        file = pd.read_csv(fname, delimiter = "\t")
        try:
            file.columns = ["rsID", "average -log10(p-value)", "latent phenotype indices"]
        except:
            pdb.set_trace()
        agg_dict = {"average -log10(p-value)": lambda x: np.mean(-np.log10(x)), 
                    "latent phenotype indices": lambda x: ", ".join(x.astype(str))}
        
        table = file.groupby('rsID').agg(agg_dict).reset_index()
        table = table.merge(pos, on = "rsID", how = "inner")
        table = table.merge(mafs, on = "rsID", how = "inner")
        colnames = ["rsID", "chromosome", "position", "minor allele frequency", 
                    "average -log10(p-value)", "latent phenotype indices"]
        table = table[colnames]
        if len(table) > 0:
            tablename = "tableS4" + letters[n] + "_" + model + "_" + E + ".txt"
            table.to_csv(tablename, sep = "\t", header = True, index = False)
            n += 1


# Figure S3 starts here
MASV_sets, feature_sets = [], []
for model in ["PCA", "NN"]:

    path_prefix = "../step7_adjust_HF_for_covariates_" + model + "/final_model_shapley_values"
    paths = [path_prefix + "/" + i for i in os.listdir(path_prefix) if "shapley_values.txt.gz" in i]
    MASVs = []
    for path in paths:
        SVs = pd.read_csv(path, compression='gzip', delimiter = "\t", low_memory = False)
        features = SVs.columns
        feature_sets.append(features[features != "eid"])
        MASVs.append(np.mean(np.abs(SVs[all_cols[all_cols != "eid"]].to_numpy()), axis = 0))

    MASV_sets.append(MASVs)

fname = "../step7_adjust_HF_for_covariates_logistic_PCA/y.txt"
features  = pd.read_csv(fname, delimiter = "\t", header = 0).columns.to_numpy()
feature_sets.append(features[features != "eid"])

logistic_PCA_path = "../step7_adjust_HF_for_covariates_logistic_PCA/logistic_SVD_output/raw_Y15.txt"
loadings = pd.read_csv(logistic_PCA_path, delimiter = "\t", header = None)
MASV_sets.append(np.abs(loadings.to_numpy()))

feature_scores = []
for i in range(3):

    MASVs, features = MASV_sets[i], feature_sets[i]
    scores = np.sum(MASVs, axis = 0)
    feature_scores.append(scores[np.argsort(features)])

# note: rankdata ranks the highest values with the highest ranks

feature_ranks = np.array([rankdata(i) for i in feature_scores]).T
rank_means = np.mean(feature_ranks, axis = 1)
sorted_mean_rank_inds = np.flip(np.argsort(rank_means))
top_features = np.sort(features)[sorted_mean_rank_inds[:62]]
top_features_rank_sets = feature_ranks[sorted_mean_rank_inds[:62]]
top_features_rank_sets = np.max(top_features_rank_sets) + 1 - top_features_rank_sets
r01 = spearmanr(top_features_rank_sets[:, 0], top_features_rank_sets[:, 1])[0]**2
r02 = spearmanr(top_features_rank_sets[:, 0], top_features_rank_sets[:, 2])[0]**2
r12 = spearmanr(top_features_rank_sets[:, 1], top_features_rank_sets[:, 2])[0]**2
figureS3a = np.concatenate([top_features.reshape(-1,1), top_features_rank_sets], axis = 1)
figureS3a = pd.DataFrame(figureS3a)
figureS3a.columns = ["ICD10 code", "PCA average importance rank", "Autoencoder average importance rank", "Logistic PCA average importance rank"]
new_codes = np.array([code.split("counts_")[-1] for code in figureS3a["ICD10 code"]])
new_codes[new_codes == "any_HF"] = "AHF"
new_codes = [code if len(code) < 4 else code[0:3] + "." + code[3:] for code in new_codes]
figureS3a["ICD10 code"] = new_codes
figureS3a.to_csv("figureS3a.txt", sep = "\t", header = True, index = False)
figureS3a_R2_vals = pd.DataFrame([[r01, r02, r12], [r01, r02, r12]])
figureS3a_R2_vals.columns = ["PCA vs NN", "PCA vs log PCA", "NN vs log PCA"]
figureS3a_R2_vals.to_csv("figureS3a_R2_vals.txt", sep = "\t", header = True, index = False)

