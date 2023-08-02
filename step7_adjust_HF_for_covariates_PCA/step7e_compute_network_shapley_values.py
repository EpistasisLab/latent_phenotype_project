import numpy as np
import pandas as pd
import os
import torch
import torch.nn as nn
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from copy import deepcopy as COPY
from tqdm import tqdm
from matplotlib import pyplot as plt
from captum.attr import IntegratedGradients
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from umap import UMAP
from time import time
import argparse
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('--index', nargs = 1, type = int, action = "store", dest = "index")
args = parser.parse_args()
index = args.index[0]
torch.set_num_threads(20)

class PCA_net(torch.nn.Module):

    # latent dims last dimension is the output shape
    def __init__(self, input_shape, output_shape):
        super(PCA_net, self).__init__()

        self.layer = nn.Linear(input_shape, output_shape)

    def forward(self, X):

        return(self.layer(X))

Y_df = pd.read_csv("y.txt", delimiter = "\t", header = 0)
ICD_cols = Y_df.columns.to_numpy()
ICD_cols = ICD_cols[ICD_cols != "eid"]
Y_np = Y_df[ICD_cols].to_numpy(dtype = int)

pcaY = (PCA(n_components = 15)).fit(Y_np)
Y_c = pcaY.transform(Y_np) 
Y_rots = (pcaY.components_).T
Y_rec = pcaY.inverse_transform(Y_c)

half_networks = [PCA_net(len(Y_np[0]), 1) for i in range(15)]
model_path = "PCA_final_model_network/final_model_PCA.model"
for i in range(15): half_networks[i].layer.weight.data = torch.tensor(Y_rots[:, i:i+1].T).float()
for i in range(15): half_networks[i].layer.bias.data = torch.zeros(1).float()

nchunks = int(len(Y_np)/5)
Y = torch.tensor(Y_np - np.mean(Y_np, axis = 0)).float()
network = half_networks[index]
IG_function = IntegratedGradients(network)
baseline = torch.mean(Y, axis = 0)
N = 3000
all_errors = torch.zeros(len(Y))
all_attributions = torch.zeros(Y.shape).float()
baseline = torch.zeros(Y.shape) + baseline
errors = torch.zeros(len(Y)).float()
attributions = torch.zeros(Y.shape).float()
boundaries = np.cumsum([0] + [int(len(Y)/nchunks)]*(nchunks - 1))
boundaries = np.array(boundaries.tolist() + [int(len(Y))])
boundaries = torch.tensor(boundaries).int()
for i in tqdm(range(nchunks)):
    out = IG_function.attribute(Y[boundaries[i]:boundaries[i + 1]], 
                                baselines = baseline[boundaries[i]:boundaries[i + 1]][0:1], 
                                method = 'riemann_trapezoid', 
                                return_convergence_delta = True, 
                                n_steps = N)
    attributions[boundaries[i]:boundaries[i + 1]] += out[0]
    errors[boundaries[i]:boundaries[i + 1]] += out[1]

fname = "final_model_shapley_values/NN_phenotype" + str(index) + "_shapley_values.txt.gz"
attribution_values = pd.DataFrame(attributions.detach().numpy())
attribution_values.columns = ICD_cols
attribution_values["eid"] = Y_df["eid"]
attribution_values.to_csv(fname, sep = "\t", header = True, index = False, compression = 'gzip')
# attribution_values2 = pd.read_csv(fname, delimiter = "\t", compression = 'gzip')

fname = "final_model_shapley_values/NN_phenotype" + str(index) + "_shapley_value_errors.txt.gz"
attribution_errors = pd.DataFrame(errors.detach().numpy().reshape(-1,1))
attribution_errors.columns = ["error"]
attribution_errors["eid"] = Y_df["eid"]
attribution_errors.to_csv(fname, sep = "\t", header = True, index = False, compression = 'gzip')
