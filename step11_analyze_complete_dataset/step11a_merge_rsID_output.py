import pandas as pd
from scipy.stats import spearmanr
import os
import pdb

base = "/home/greggj/pleiotropy_and_GxE/step10_get_significant_SNPs_"
# prefixes = [base + i for i in ["NN", "PCA", "logistic_PCA"]]
prefix_ends = ["NN", "PCA", "logistic_PCA"]
prefixes = [base + i for i in prefix_ends]
suffixes = ["GxAlcohol", "GxExercise", "GxSmoking", "GxGender", "main"]
suffixes = ["rsIDs_" + i + "_effects.txt" for i in suffixes]

bim_suffixes = ["alcohol", "exercise", "smoking", "gender"]
bim_suffixes = ["significant_SNPs_plink_files_" + i for i in bim_suffixes]
bim_files = []
table1b = []
for prefix in prefixes:
    table1b_path = prefix + "/step10e_p_val_analysis.txt"
    table1b.append(pd.read_csv(table1b_path, delimiter = "\t")) 

    paths = [prefix + "/" + suffix for suffix in bim_suffixes]
    all_paths = []
    for path in paths: all_paths += [path + "/" + i for i in os.listdir(path)]
    all_paths = [i for i in all_paths if ".bim" in i]
    info = [pd.read_csv(path, delimiter = "\t", header = None) for path in all_paths]
    bim_files.append(pd.concat(info)[[1, 0, 3]])

table1b = pd.concat(table1b)
table1b["model"] = prefix_ends
ordered_names = ["model", "TRACE vs. linear logp Mean", "p-value", "conditional logp Mean", "conditional p-value"]
new_names = ["latent phenotype model", "TRACE vs. linear E[logp]", "p-value", "conditional E[logp]", "conditional p-value"]
table1b = table1b[ordered_names]
table1b.columns = new_names
table1b.to_csv("table1b.txt", sep = "\t", header = True, index = False)

bim_files = pd.concat(bim_files).drop_duplicates(1)
bim_files.columns = ["rsID", "chr", "pos"]

QTL_suffixes = ["alcohol", "exercise", "smoking", "gender"]
QTL_suffixes = ["hits_QTL_" + i for i in QTL_suffixes]
QTL_files = []
for prefix, end in zip(prefixes, prefix_ends):
    paths = [prefix + "/" + i for i in QTL_suffixes]
    files = []
    for path in paths:
        env_factor = path.split("_")[-1]
        all_paths = [path + "/" + i for i in os.listdir(path)]
        file = pd.concat([pd.read_csv(path, delimiter = "\t") for path in all_paths])
        file["type"] = end
        file["env_factor"] = env_factor
        files.append(file)
    QTL_files.append(pd.concat(files))

QTL_files = pd.concat(QTL_files)[["rsID", "pheno_index", "p_main", "p_null2", "type", "env_factor"]]
QTL_files["p"] = (QTL_files["p_main"]/QTL_files["p_null2"]).to_numpy()
QTL_files = QTL_files[["rsID", "pheno_index", "p", "type", "env_factor"]]

frq_base = "/home/greggj/pleiotropy_and_GxE/step9_regress_phenotypes_against_SNPs_NN/genotype_metadata/"
chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "XY", "Y", "MT" ]
frq_paths = [frq_base + "genotype_metadata_chr" + chr + ".frq" for chr in chromosomes]
frq_files = pd.concat([pd.read_csv(name, delim_whitespace = True, usecols = ["SNP", "MAF"]) for name in frq_paths])
frq_files.columns = ["rsID", "MAF"]

test_efficacy_dfs = []
factors = ["alcohol", "exercise", "smoking", "gender", "main"]
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

pvals = pd.concat(test_efficacy_dfs)[["p", "pEDGE2"]]
spearmanr(pvals["p"], pvals["pEDGE2"])
