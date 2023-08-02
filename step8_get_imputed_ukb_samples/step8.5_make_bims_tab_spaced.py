import pandas as pd
import pdb

prefix = "filtered_output/UKB_samples_chr"
chr_vals = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chr_vals += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
paths = [prefix + i + ".bim" for i in chr_vals]

for path in paths:
    file = pd.read_csv(path, delim_whitespace = True, header = None)
    file.to_csv(path, sep = "\t", header = False, index = False)