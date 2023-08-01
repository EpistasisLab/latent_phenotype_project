import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.stats import mannwhitneyu as MWU
from functools import reduce
from tqdm import tqdm
import os
import pdb

base = "/home/greggj/pleiotropy_and_GxE/step10_get_significant_SNPs_"
# prefixes = [base + i for i in ["NN", "PCA", "logistic_PCA"]]
prefix_ends = ["NN", "PCA", "logistic_PCA"]
prefixes = [base + i for i in prefix_ends]
suffixes = ["GxAlcohol", "GxExercise", "GxSmoking", "GxGender", "main"]
suffixes = ["rsIDs_" + i + "_effects.txt" for i in suffixes]
files = []
for prefix in prefixes:
    files_model = []
    for suffix in suffixes:
        path = prefix + "/" + suffix
        file = pd.read_csv(path, delimiter = "\t")
        files_model.append(pd.read_csv(path, delimiter = "\t"))
    files.append(files_model)

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "XY", "X"]
step9_prefix = "../step9_regress_phenotypes_against_SNPs_"
step9_prefixes = [step9_prefix + end + "/QTL_output_" for end in prefix_ends]
for prefix, p_start in zip(step9_prefixes, ["NN", "PCA", "logistic_PCA"]):
    for k, env in enumerate(["alcohol", "exercise", "smoking", "gender", "gender"]):
        if k == 4: dir = p_start + "_" + "main"
        else: dir = p_start + "_" + env
        if not os.path.isdir(dir):
            os.mkdir(dir)
        if not os.path.exists(dir + "/" + "pvals.txt"):
            min_p_val_sets = []
            for chr in tqdm(chromosomes):
                if p_start == "NN":
                    Np = 15
                else:
                    Np = 16
                paths = [prefix + env + "/" + "QTL_effects_chr" + chr + "_P" + str(ind) + ".txt" for ind in range(Np)]
                if k == 4:
                    data = [pd.read_csv(path, delimiter = "\t", usecols = [0, 3]) for path in paths]
                    for i in range(Np): data[i].columns = [0, i + 1]
                else:
                    data = [pd.read_csv(path, delimiter = "\t", usecols = [0, 1, 3]) for path in paths]
                    for i in range(Np): data[i]["proxy"] = data[i]["p_main"]/data[i]["p_null2"] 
                    for i in range(Np): data[i] = data[i][["rsID", "proxy"]] 
                    for i in range(Np): data[i].columns = [0, i + 1]
                def ij(df1, df2): return(df1.merge(df2, on = 0, how = "inner"))
                all_p_vals = reduce(ij, data)[np.arange(1, Np + 1)].to_numpy()
                min_p_val_sets.append(np.min(all_p_vals, axis = 1))           
            min_p_val_sets = np.concatenate(min_p_val_sets)
            min_p_val_sets = min_p_val_sets[min_p_val_sets != 0]
            df = pd.DataFrame(-np.log10(min_p_val_sets))
            df.to_csv(dir + "/" + "pvals.txt", sep = "\t", header = True, index = False)


stats_main, stats_GxSmoking = [], []
for pref, pref_end, files_model in zip(step9_prefixes, prefix_ends, files):
    # note that main effects are documented in all 4 paths.
    # main effects correspond to the second "gender" in step9_paths
    # and will be handled accordingly
    step9_paths = [pref + i for i in ["alcohol", "exercise", "smoking", "gender", "gender"]]
    for main_path, file, effect in zip(step9_paths, files_model, suffixes):
        file.columns = ["rsID"] + [name + "_" + pref_end + "_step10" for name in file.columns if name != "rsID"]
        for model in prefix_ends:
            if model == pref_end:
                continue
            if (effect != "rsIDs_main_effects.txt") and (effect != "rsIDs_GxSmoking_effects.txt"):
                continue
            model_path = main_path.replace(pref_end, model) 
            # phenotype 15 for the NN, the residual phenotype, contains less than 0.01% of the variance. 
            if model == "NN": paths = [model_path + "/" + i for i in os.listdir(model_path) if "P15" not in i]      
            else: paths = [model_path + "/" + i for i in os.listdir(model_path)]
            step9_files = []
            if effect == "rsIDs_main_effects.txt":
                for path in tqdm(paths):
                    pheno_index = int(path.split("_P")[-1].split(".")[0])
                    all_SNPs = pd.read_csv(path, delimiter = "\t", usecols = [0, 3], low_memory = False)
                    proxy_p_name = 'p_null2_' + model + "_step9_proxy"
                    main_p_name = 'p_null2_' + pref_end + "_step10"
                    all_SNPs.columns = ["rsID", proxy_p_name]
                    all_SNPs['pheno_index_' + model + '_step9'] = pheno_index
                    all_SNPs = all_SNPs.merge(file, on = "rsID", how = "inner")
                    step9_files.append(all_SNPs)
            else: 
                for path in tqdm(paths):
                    pheno_index = int(path.split("_P")[-1].split(".")[0])
                    all_SNPs = pd.read_csv(path, delimiter = "\t", usecols = [0, 1, 3], low_memory = False)
                    all_SNPs['pEDGE2_' + model + "_step9_proxy"] = all_SNPs["p_main"]/all_SNPs["p_null2"]
                    proxy_p_name = 'pEDGE2_' + model + "_step9_proxy"
                    main_p_name = 'pEDGE2_' + pref_end + "_step10"
                    all_SNPs = all_SNPs[["rsID", proxy_p_name]]
                    all_SNPs['pheno_index_' + model + '_step9'] = pheno_index
                    all_SNPs = all_SNPs.merge(file, on = "rsID", how = "inner")
                    step9_files.append(all_SNPs)
            step9_files = pd.concat(step9_files)
            closest_matches = []
            for rsID in np.unique(file["rsID"].to_numpy()):
                num_effects = len(file.loc[file["rsID"] == rsID, :])
                subset = step9_files.loc[step9_files["rsID"] == rsID, :]
                subset = subset.reset_index()
                del subset["index"]
                inds = np.argsort(subset[proxy_p_name].to_numpy())
                closest_matches.append(subset.loc[inds[:num_effects], :])
            closest_matches = pd.concat(closest_matches)
            alt_p_vals = -np.log10(closest_matches[proxy_p_name].to_numpy())
            real_p_vals = -np.log10(closest_matches[main_p_name].to_numpy())
            if "main_effects" in effect: env_factor = "main"
            else: env_factor = main_path.split("_")[-1]
            data_name = "realp_" + pref_end + "_altp_" + model + "_" + env_factor
            null_p_fname = model + "_" + env_factor + "/" + "pvals.txt"
            null_p_vals = pd.read_csv(null_p_fname, delimiter = "\t").to_numpy().reshape(-1)
            p_diff = MWU(alt_p_vals, null_p_vals, alternative = "greater")[1]
            output = [data_name, np.mean(real_p_vals), np.mean(alt_p_vals), np.mean(null_p_vals), p_diff]
            if effect == "rsIDs_main_effects.txt": stats_main.append(output)
            elif effect == "rsIDs_GxSmoking_effects.txt": stats_GxSmoking.append(output)

df_main = pd.DataFrame(stats_main)
df_main = df_main.sort_values(0, ascending = False)
df_main [5] = ["logistic PCA", "logistic PCA", "PCA", "PCA", "autoencoder", "autoencoder"]
df_main [6] = ["PCA", "autoencoder", "logistic PCA", "autoencoder", "logistic PCA", "PCA"]
df_main = df_main[[5, 6, 1, 2, 3, 4]] 
df_main.columns = ["main model", "alternate model",  "main E[-log10(p)]", "alt E[-log10(p)]", "null E[-log10(p)]", "p_diff alt minus null"]
df_main.to_csv("table1c.txt", sep = "\t", header = True, index = False)

df_smoking = pd.DataFrame(stats_GxSmoking)
df_smoking = df_smoking.sort_values(0, ascending = False)
df_smoking [5] = ["logistic PCA", "logistic PCA", "PCA", "PCA", "autoencoder", "autoencoder"]
df_smoking [6] = ["PCA", "autoencoder", "logistic PCA", "autoencoder", "logistic PCA", "PCA"]
df_smoking = df_smoking[[5, 6, 1, 2, 3, 4]] 
df_smoking.columns = ["main model", "alternate model",  "main E[-log10(p)]", "alt E[-log10(p)]", "null E[-log10(p)]", "p_diff alt minus null"]
df_smoking.to_csv("table1d.txt", sep = "\t", header = True, index = False)
