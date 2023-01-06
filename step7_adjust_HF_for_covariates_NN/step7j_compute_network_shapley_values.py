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

X_df = pd.read_csv("X.txt", delimiter = "\t", header = 0)
X_cols = np.array(X_df.columns)
X_cols[X_cols == "22001-0.0"] = "is_male"
X_df.columns = X_cols
X = X_df[X_cols[X_cols != "eid"]].to_numpy() 
X_cols = X_cols[X_cols != "eid"]
is_PC = np.array([col[0:2] == "PC" for col in X_cols])
PCs = X[:, is_PC]

prefix = "../step1_get_phenotypes_simple/"
death_ICD_data_all = pd.read_csv(prefix + "death_ICD_data.txt", delimiter = "\t", low_memory = False)
filtered_eids = pd.read_csv("y.txt", delimiter = "\t", usecols = ["eid"])
death_ICD_data_all = death_ICD_data_all.merge(filtered_eids, on = "eid", how = "inner")
data_cols = death_ICD_data_all.columns.to_numpy()[death_ICD_data_all.columns.to_numpy() != "eid"]
death_ICD_data = death_ICD_data_all[data_cols]

ICD_seqs = [row[pd.isnull(row) == False].tolist() for row in death_ICD_data.to_numpy()]
ICD_seqs = [row if len(row) > 0 else [""] for row in ICD_seqs]
ICD_codes = death_ICD_data.to_numpy().reshape(-1)
ICD_codes = ICD_codes[pd.isnull(ICD_codes) == False]
ICD_codes = np.unique(ICD_codes)

ICD_codes2 = np.array(ICD_codes).reshape(-1, 1)
Y = np.zeros((len(ICD_seqs), len(ICD_codes)), dtype = bool)
for i, seq in enumerate(ICD_seqs): Y[i, :] = np.any(np.array(seq).reshape(-1, 1) == np.array(ICD_codes), axis = 0)
out_dim = len(Y[0])

main_network = simpleAE(len(Y[0]), [1000, 500, 15], 0)
half_networks = [half_AE(len(Y[0]), [1000, 500, 1], 0) for i in range(15)]
main_network.eval()
for half_net in half_networks: half_net.eval()
model_path = "normalAE_final_model_fold_networks/final_model_normalAE_0.3dop_0rs_5ws_1000d1_500d2.model"
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

pheno_fname = "step7i_0.9628510715767036r_output_latent_phenotypes.txt"
raw_phenotypes = pd.read_csv(pheno_fname, delimiter = "\t").to_numpy()
if index == 0:
    for i in range(15): 
        pheno = raw_phenotypes[:, i]
        vals, counts = np.unique(pheno, return_counts = True)
        spikes = vals[np.flip(np.argsort(counts))][0:2]
        pheno = pheno[np.isin(pheno, spikes) == False] 
        plt.hist(pheno, 200)
        plt.savefig("step7j_pheno" + str(i) + ".png")
        plt.clf()

nchunks = int(len(Y)/3)
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

fname1 = "final_model_shapley_values/step7j_network" + str(index) + "_shapley_values.txt"
fname2 = "final_model_shapley_values/step7j_network" + str(index) + "_shapley_value_errors.txt"
pd.DataFrame(attributions.detach().numpy()).to_csv(fname1, sep = "\t", header = False, index = False)
pd.DataFrame(errors.detach().numpy().reshape(-1,1)).to_csv(fname2, sep = "\t", header = False, index = False)

'''

#TODO: 
#1) create partial neural network with Y input and Z output
#2) analyze Y->Z shapley values
#3) reconstruct HF phenotype and get residuals
#4) correct for genetic PCs

HF_residuals = Y[:, 4] - Y_rec[:, 4]
Y_c_res = np.concatenate([Y_c, HF_residuals.reshape(-1, 1)], axis = 1)
PCs2 = np.concatenate((PCs, np.ones((len(PCs), 1))), axis = 1)
PCsT_PCs = np.matmul(PCs2.T, PCs2)
correction = np.eye(len(PCsT_PCs))*np.min(PCsT_PCs)*1E-6
PCsT_PCs_inv = np.linalg.inv(PCsT_PCs + correction)
coef = np.matmul(PCsT_PCs_inv, PCs2.T)
weights = np.matmul(coef, Y_c_res)
Y_est = np.matmul(PCs2, weights)
phenotypes = Y_c_res - Y_est

# put fake phenotype generation code here if replacing

phenotypes = pd.DataFrame(phenotypes)
phenotypes["eid"] = X_df["eid"].astype(int)
phenotypes = phenotypes[["eid"] + phenotypes.columns[phenotypes.columns != "eid"].tolist()]
phenotypes.to_csv("logistic_SVD_output/phenotypes" + str(k) + ".txt", sep = "\t", header = True, index = False)
'''

