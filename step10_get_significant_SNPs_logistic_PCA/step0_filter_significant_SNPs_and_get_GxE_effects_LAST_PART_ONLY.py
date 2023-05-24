import numpy as np
import pandas as pd
import argparse
import os
import psutil
from functools import reduce
from tqdm import tqdm
from copy import deepcopy as COPY
from scipy.stats import linregress
import statsmodels.api as sm
from scipy.stats import chi2
from bed_reader import open_bed
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
import pdb
import matplotlib as mpl
from scipy.stats import norm
mpl.rcParams['agg.path.chunksize'] = 10000
from statsmodels.stats.outliers_influence import summary_table

parser = argparse.ArgumentParser()
parser.add_argument('--chr', nargs = 1, type = int, action = "store", dest = "chr")
parser.add_argument('--name', nargs = 1, type = str, action = "store", dest = "name")
args = parser.parse_args()
chr_index = (args.chr)[0] - 1
name = (args.name)[0]

if not os.path.exists("hits_QTL" + name):
    os.mkdir("hits_QTL" + name)
if not os.path.exists("hits_GxE_p_vals_getters" + name):
    os.mkdir("hits_GxE_p_vals_getters" + name)

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]

path = "hits_QTL" + name + "/QTL_hits_chr" + chromosomes[chr_index] + ".txt"
best_rsIDs3 = pd.read_csv(path, delimiter = "\t")

for i in range(len(best_rsIDs3)):
    info_i = best_rsIDs3.loc[i, ["rsID", "pheno_index", "chr"]].to_numpy(dtype = str)
    name_parts = COPY(info_i)
    name_parts[0] = "_".join(name_parts[0].split(":"))
    file = open("hits_GxE_p_vals_getters" + name + "/" + "_".join(name_parts) + ".sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#BSUB -J GxE_p_getter\n")
    file.write("#BSUB -o hits_GxE_p_vals_getters" + name + "/" + "_".join(name_parts) + ".out\n")
    file.write("#BSUB -e hits_GxE_p_vals_getters" + name + "/" + "_".join(name_parts) + ".err\n")
    file.write("source activate torch_env2\n\n")
    file.write("python step0_compute_GxE_p_values.py ")
    file.write("--rsID " + info_i[0] + " --pheno_index " + str(info_i[1]) + " --chr " + info_i[2] + " --name " + name)
    file.close()

pdb.set_trace()
print(1)