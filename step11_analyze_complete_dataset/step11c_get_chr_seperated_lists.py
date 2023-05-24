import pandas as pd
import os
import pdb

if not os.path.isdir("snps"):
    os.mkdir("snps")
if not os.path.isdir("snps_output"):
    os.mkdir("snps_output")

for effect in ["main_effects", "GxSmoking_effects", "GxAlcohol_effects", "GxGender_effects"]:
    file = pd.read_csv("rsIDs_" + effect + ".txt", delimiter = "\t").astype(str)
    chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
    chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "XY"]
    subsets = [file[file["chr"] == ch] for ch in chromosomes if len(file[file["chr"] == ch]) > 0]
    for subset in subsets:
        snps = subset[["rsID"]].drop_duplicates().to_numpy().reshape(-1)
        for snp in snps:
            path = "snps/" + effect + "_" + snp + ".txt"
            path = path.replace(":", "!")
            snps = pd.DataFrame([snp])
            snps.to_csv(path, sep = "\t", header = False, index = False)
