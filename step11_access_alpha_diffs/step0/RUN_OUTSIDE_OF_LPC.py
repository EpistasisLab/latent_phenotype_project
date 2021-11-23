import numpy as np 
import pandas as pd
import tabula
import pdb
import requests
import sys
 
server = "https://rest.ensembl.org"
chromosomes = list(np.arange(1, 23).astype(str)) + ["X", "Y", "MT"]
ext_vec = ["/map/human/GRCh37/" + str(i) + "/GRCh38?" for i in chromosomes]
headers = { "Content-Type" : "application/json"}
r_vec = [requests.get(server+ext, headers=headers) for ext in ext_vec]
decoded = [rr.json()["mappings"] for rr in r_vec]
# "mapped" is GRCh38, while "original" is GRCh37
for i, d in enumerate(decoded):
    mat = np.zeros((len(d), 5))
    for j, r in enumerate(d):
        mat[j, 0] = r["original"]["start"]
        mat[j, 1] = r["original"]["end"]
        mat[j, 2] = r["mapped"]["start"]
        mat[j, 3] = r["mapped"]["end"]
        mat[j, 4] = r["mapped"]["start"] - r["original"]["start"]
    if np.any(mat[:, 1] - mat[:, 0] != mat[:, 3] - mat[:, 2]):
        print("exiting: the mapping format is wrong or unexpected")
        exit()
    df = pd.DataFrame(mat)
    df.columns = ["hg37_start", "hg37_end", "hg38_start", "hg38_end", "shift"]
    fname = "chr" + chromosomes[i] + "_mapper.txt"
    df.to_csv(fname, sep = "\t", header = False, index = False)

