import pandas as pd
from scipy.stats import spearmanr
import os
import pdb

base = "/home/greggj/pleiotropy_and_GxE/step10_get_significant_SNPs_"
# prefixes = [base + i for i in ["NN", "PCA", "logistic_PCA"]]
prefix_ends = ["NN", "PCA", "logistic_PCA"]
prefixes = [base + i for i in prefix_ends]
suffixes = ["main"] #["GxAlcohol", "GxExercise", "GxSmoking", "GxGender", "main"]
suffixes = ["rsIDs_" + i + "_effects.txt" for i in suffixes]

bim_suffixes = ["alcohol", "exercise", "smoking", "gender"]
bim_suffixes = ["significant_SNPs_plink_files_" + i for i in bim_suffixes]
bim_files = []
for prefix in prefixes:
    paths = [prefix + "/" + suffix for suffix in bim_suffixes]
    all_paths = []
    for path in paths: all_paths += [path + "/" + i for i in os.listdir(path)]
    all_paths = [i for i in all_paths if ".bim" in i]
    info = [pd.read_csv(path, delimiter = "\t", header = None) for path in all_paths]
    bim_files.append(pd.concat(info)[[1, 0, 3]])

bim_files = pd.concat(bim_files).drop_duplicates(1)
bim_files.columns = ["rsID", "chr", "pos"]

frq_base = "/home/greggj/pleiotropy_and_GxE/step9_regress_phenotypes_against_SNPs_PCA/genotype_metadata/"
chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "XY", "Y", "MT" ]
frq_paths = [frq_base + "genotype_metadata_chr" + chr + ".frq" for chr in chromosomes]
frq_files = pd.concat([pd.read_csv(name, delim_whitespace = True, usecols = ["SNP", "MAF"]) for name in frq_paths])
frq_files.columns = ["rsID", "MAF"]

test_efficacy_dfs = []
factors = ["main"] # ["alcohol", "exercise", "smoking", "gender", "main"]
for suffix, env_factor in zip(suffixes, factors):
    merged_rsID_set = []
    for prefix, end in zip(prefixes, prefix_ends):

        path = prefix + "/" + suffix
        df = pd.read_csv(path, delimiter = "\t")

        folder = "all_sig_rsIDs_" + prefix.split("SNPs_")[-1]
        if not os.path.isdir(folder): os.mkdir(folder)
        path2 = folder + "/" + suffix
        df.to_csv(path2, sep = "\t", header = True, index = False)

        df["type"] = end
        df["env_factor"] = env_factor
        merged_rsID_set.append(df)


    if suffix != "rsIDs_main_effects.txt":
        test_efficacy_df = pd.concat(merged_rsID_set).merge(QTL_files, on = ["rsID", "pheno_index", "type", "env_factor"], how = "inner")
        test_efficacy_dfs.append(test_efficacy_df)
        df = pd.concat(merged_rsID_set)[["rsID", "pEDGE2"]].drop_duplicates()
    else:
        df = pd.concat(merged_rsID_set)[["rsID", "p_null2"]].drop_duplicates()

    df.columns = ["rsID", "P-value"]
    SNP_MAFs = df.merge(frq_files, on = "rsID", how = "inner")
    SNP_MAFs.to_csv("SNP_MAFs_" + suffix, sep = "\t", header = True, index = False)
    df.to_csv(suffix[:-4] + "_pvals.txt", sep = "\t", header = True, index = False)
    old_len = len(df)
    df = df.merge(bim_files, on = "rsID", how = "inner")
    new_len = len(df)
    if old_len != new_len:
        pdb.set_trace()
    df[["rsID", "chr", "pos"]].to_csv(suffix, sep = "\t", header = True, index = False)