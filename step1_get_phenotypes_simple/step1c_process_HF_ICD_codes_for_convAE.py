import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
from matplotlib import pyplot as plt
from copy import deepcopy as COPY
from functools import reduce
from tqdm import tqdm
import time
import pdb
import os

fields = pd.read_csv("phenotype_info.txt", delimiter = "\t", header = 0,  dtype = str)

is_white = np.isin(fields["21000-0.0"], ["1", "1001", "1002", "1003"])
self_declared_illness_cols = fields.columns[np.array(["20002" in col for col in fields.columns])]
has_HCM1 = np.any(fields[self_declared_illness_cols] == "1588" , axis = 1)
has_HCM2 = np.any(np.isin(fields, ["I421", "I422"]), axis = 1)
has_no_HCM = np.logical_or(has_HCM1, has_HCM2) == False
fields = fields.loc[np.logical_and(is_white, has_no_HCM), :]
fields[["eid", "eid"]].to_csv("../step2_get_UKB_samples/eids.tab", sep = "\t", header = False, index = False)
colnames = fields.columns
colfields = np.array([name.split("-")[0] for name in colnames])

binary_ICD_data = pd.read_csv("Y.txt", delimiter = "\t", header = 0)
normal_ICD_data = fields[colnames[colfields == "41270"]]
normal_ICD_mat = normal_ICD_data.to_numpy()

death1_data = fields[colnames[colfields == "40001"]].to_numpy()
z1 = COPY(death1_data)
for i in tqdm(range(len(death1_data))): 
    if death1_data[i, 0] in normal_ICD_mat[i]: 
        death1_data[i, 0] = np.nan
    
death2_data = fields[colnames[colfields == "40002"]].to_numpy()
z2 = COPY(death2_data)
for i in tqdm(range(len(death2_data))): 
    already_happened = np.isin(death2_data[i], normal_ICD_mat[i])
    death2_data[i, already_happened] = np.nan

death_ICD_data = pd.DataFrame(np.concatenate([normal_ICD_mat, death1_data, death2_data], axis = 1))
my_ICD_codes = [i.split("_")[1] for i in binary_ICD_data.columns[1:]]
my_ICD_codes += ["I110", "I130", "I132", "I255", "I420", "I425", "I426", "I427", "I428", "I429", "I431", "I500", "I501", "I509"]
in_my_ICD_codes = death_ICD_data.isin(my_ICD_codes)
death_ICD_data[in_my_ICD_codes == False] = np.nan
death_ICD_data["eid"] = fields["eid"].to_numpy()
death_ICD_data.to_csv("death_ICD_data.txt", sep = "\t", header = True, index = False)

'''
ICD_seqs = [row[pd.isnull(row) == False] for row in death_ICD_data.to_numpy()]

seq_lengths = [len(row) for row in ICD_seqs]
L_max = np.max(seq_lengths)

my_ICD_codes = np.array(my_ICD_codes).reshape(-1, 1)
ICD_ind_seqs = [np.where(np.any(my_ICD_codes ==  seq, axis = 1))[0] for seq in ICD_seqs]
one_hot_seqs = [np.ones((L_max, len(my_ICD_codes)))*np.nan for l in seq_lengths]
for i in tqdm(range(len(one_hot_seqs))): 
    for k in range(seq_lengths[i]):       
        m = ICD_ind_seqs[i][k]
        one_hot_seqs[i][k][m] = 1

Y_convAE = pd.DataFrame(np.concatenate(one_hot_seqs, axis = 0))
inds = np.arange(len(one_hot_seqs[0][0])).tolist()
del one_hot_seqs
eids = np.repeat(fields["eid"].to_numpy(), L_max).astype(int)
Y_convAE.columns = inds
Y_convAE["eid"] = eids
Y_convAE.to_csv("Y_convAE.txt", sep = "\t", header = True, index = False)
'''
