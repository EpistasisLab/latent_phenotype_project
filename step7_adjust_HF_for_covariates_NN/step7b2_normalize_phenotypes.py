import numpy as np
import pandas as pd
import os
import pdb 
import argparse
from tqdm import tqdm
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import yeojohnson as yj
from copy import deepcopy as COPY

Y = pd.read_csv("y.txt", delimiter = "\t")
ICD_coders = np.sum(Y.to_numpy(dtype = int)[:, 1:], axis = 1) > 1
P = pd.read_csv("phenotypes_NN_R0.9620560038092117.txt", delimiter = "\t")
P2 = COPY(P)
for i in np.arange(16).astype(str):
    pheno = (P[i].to_numpy() - np.mean(P[i].to_numpy()))/np.std(P[i].to_numpy())
    pheno2 = yj(pheno)[0]
    P2.loc[:, i] = pheno2
P2.to_csv("phenotypes_for_step9.txt", sep = "\t", header = True, index = False)
