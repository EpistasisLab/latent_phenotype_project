import numpy as np
import pandas as pd
import statsmodels.api as sm
import os
import pdb
from copy import deepcopy as COPY
from bed_reader import open_bed
from functools import reduce
from scipy.stats import linregress

paths = ["hits_QTL_gender"]
GxE_paths = ["hits_GxE_p_vals"]
for path, GxE_path in zip(paths, GxE_paths):
    pdb.set_trace()
    col4 = [name[:-4].split("_")[-1] for name in os.listdir(GxE_path) if ".txt" in name]
    col3 = [name[:-4].split("_")[-2] for name in os.listdir(GxE_path) if ".txt" in name]
    col2 = [name[:-4].split("_")[-3] for name in os.listdir(GxE_path) if ".txt" in name]
    col1 = ['_'.join(name[:-4].split("_")[:-3]) for name in os.listdir() if ".txt" in name]
    data = pd.DataFrame(np.array([col1, col2, col3, col4]).T)
    data[3] = data.loc[:, 3].astype(float)
    data = data.sort_values(by = 3, ascending = True)
    data.to_csv("check1.txt", sep = "\t", header = False, index = False)

