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

class simpleAE(torch.nn.Module):

    # latent dims last dimension is the output shape
    def __init__(self, output_shape, latent_dims, dop):
        super(simpleAE, self).__init__()

        dims_in = [output_shape] + latent_dims
        dims_out = np.flip(dims_in)
        self.in_layers = nn.Sequential()
        self.out_layers = nn.Sequential()
        
        for i in range(len(dims_in) - 1):
             self.in_layers.add_module("in_relu" + str(i), torch.nn.LeakyReLU(negative_slope = 0.05))
             if i > 0: self.in_layers.add_module("in_dropout" + str(i), nn.Dropout(p = dop))
             self.in_layers.add_module("in_linear" + str(i + 1), nn.Linear(dims_in[i], dims_in[i + 1]))
        for i in range(len(dims_in) - 1):
             layer_i = nn.Linear(dims_out[i], dims_out[i + 1])
             self.out_layers.add_module("out_relu" + str(i), torch.nn.LeakyReLU(negative_slope = 0.05))
             if i > 0: self.out_layers.add_module("out_dropout" + str(i), nn.Dropout(p = dop))
             self.out_layers.add_module("out_linear" + str(i + 1), layer_i)

        self.sig = nn.Sigmoid()

    def forward(self, X):

        Z = self.in_layers(X)
        Y_pred = self.sig(self.out_layers(Z))
        # Y_pred = self.out_layers(Z)
        return(Z, Y_pred)

class half_AE(torch.nn.Module):

    # latent dims last dimension is the output shape
    def __init__(self, output_shape, latent_dims, dop):
        super(half_AE, self).__init__()

        dims_in = [output_shape] + latent_dims
        dims_out = np.flip(dims_in)
        self.in_layers = nn.Sequential()
        
        for i in range(len(dims_in) - 1):
             self.in_layers.add_module("in_relu" + str(i), torch.nn.LeakyReLU(negative_slope = 0.05))
             if i > 0: self.in_layers.add_module("in_dropout" + str(i), nn.Dropout(p = dop))
             self.in_layers.add_module("in_linear" + str(i + 1), nn.Linear(dims_in[i], dims_in[i + 1]))

    def forward(self, X):

        Z = self.in_layers(X)
        return(Z)

Y_df = pd.read_csv("y.txt", delimiter = "\t", header = 0)
ICD_cols = Y_df.columns.to_numpy()
ICD_cols = ICD_cols[ICD_cols != "eid"]
Y = Y_df[ICD_cols].to_numpy(dtype = int)

main_network = simpleAE(len(Y[0]), [1000, 500, 15], 0)
half_networks = [half_AE(len(Y[0]), [1000, 500, 1], 0) for i in range(15)]
main_network.eval()
for half_net in half_networks: half_net.eval()
model_path = "AE_final_model_network/final_model_AE_0.3dop_0rs_1000d1_500d2.model"
main_network.load_state_dict(torch.load(model_path))
param_list = [mat.data for mat in main_network.parameters()][:6]
for i in range(15): half_networks[i].in_layers[1].weight.data = param_list[0]
for i in range(15): half_networks[i].in_layers[1].bias.data = param_list[1]
for i in range(15): half_networks[i].in_layers[4].weight.data = param_list[2]
for i in range(15): half_networks[i].in_layers[4].bias.data = param_list[3]
for i in range(15): half_networks[i].in_layers[7].weight.data = param_list[4][i:i+1]
for i in range(15): half_networks[i].in_layers[7].bias.data = param_list[5][i:i+1]

'''
correct_comps = []
for i in range(15): correct_comps.append(torch.all(torch.abs(main_network(Y)[0][:, i:i+1] - half_networks[i](Y)) < 0.0001))
torch.all(torch.tensor(correct_comps))
'''

nchunks = int(len(Y)/5)
Y = torch.tensor(Y).float()
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

fname = "final_model_shapley_values/NN_phenotype" + str(index) + "_shapley_value_errors.txt.gz"
attribution_errors = pd.DataFrame(errors.detach().numpy().reshape(-1,1))
attribution_errors.columns = ["error"]
attribution_errors["eid"] = Y_df["eid"]
attribution_errors.to_csv(fname, sep = "\t", header = True, index = False, compression = 'gzip')
