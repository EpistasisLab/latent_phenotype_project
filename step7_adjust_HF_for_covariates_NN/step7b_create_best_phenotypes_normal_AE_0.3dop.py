import numpy as np
import pandas as pd
import os
import torch
import torch.nn as nn
import argparse
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from sklearn.model_selection import KFold
from copy import deepcopy as COPY
from tqdm import tqdm
from matplotlib import pyplot as plt
from gensim.test.utils import common_texts
from gensim.models import Word2Vec
from gensim.models.callbacks import CallbackAny2Vec
from sklearn.cluster import KMeans
from umap import UMAP
from time import time
from itertools import chain
import pdb
import sys
import gc
import psutil

parser = argparse.ArgumentParser()
parser.add_argument('--d1', nargs = 1, type = int, action = "store", dest = "dim1")
parser.add_argument('--rs', nargs = 1, type = int, action = "store", dest = "rs")
parser.add_argument('--dop', nargs = 1, type = float, action = "store", dest = "dop")

args = parser.parse_args()
rs = args.rs[0]
dim1 = args.dim1[0]
dim2 = int(dim1/2)
dop = args.dop[0]

if not (os.path.exists("X.txt") and os.path.exists("y.txt")):
    X = pd.read_csv("../step1_get_phenotypes_simple/X.txt", delimiter = "\t", header = 0)
    y = pd.read_csv("../step1_get_phenotypes_simple/y.txt", delimiter = "\t", header = 0)
    Y = pd.read_csv("../step1_get_phenotypes_simple/Y.txt", delimiter = "\t", header = 0)
    PCs = pd.read_csv("../step6_PCA/UKB_samples_unrelated_pruned.evec", skiprows = 1, delim_whitespace = True, header = None)
    PCs.columns = ["eid"] + ["PC" + str(i) for i in range(1, len(PCs.columns))]
    del PCs[PCs.columns[-1]]
    unrelated_eids = pd.read_csv("../step4_remove_relatives/UKB_samples_unrelated.fam", usecols = [0], delimiter = " ", header = None)
    unrelated_eids.columns = ["eid"]

    X = X.merge(unrelated_eids, on = "eid", how = "inner")
    X = X.merge(PCs, on = "eid", how = "inner")
    X.to_csv("X.txt", sep = "\t", header = True, index = False)

    y = y.merge(unrelated_eids, on = "eid", how = "inner")
    y = y.merge(Y, on = "eid", how = "inner")
    y.to_csv("y.txt", sep = "\t", header = True, index = False)

X_df = pd.read_csv("X.txt", delimiter = "\t", header = 0)
X_cols = np.array(X_df.columns)
X_cols = X_cols[X_cols != "eid"]
X = X_df[X_cols].to_numpy() 
is_PC = np.array([col[0:2] == "PC" for col in X_cols])
PCs = X[:, is_PC]

# prints ram usage (GB)
process = psutil.Process(os.getpid())
print("baseline RAM (GB):" + str(process.memory_info().rss/1E9))
torch.set_num_threads(20)
torch.manual_seed(0)
np.random.seed(0)

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

def get_batches(X_train, X_test, Y_train, Y_test, batch_size):

    sample_size_train = len(X_train)
    num_divisible_train = (sample_size_train - sample_size_train%batch_size)
    num_batches = int(num_divisible_train/batch_size)
    batch_indices_train = np.arange(sample_size_train)
    np.random.shuffle(batch_indices_train)
    batch_indices_train = torch.tensor(batch_indices_train).long()
    batch_indices_train = batch_indices_train[:num_divisible_train]

    sample_size_test = len(X_test)
    batch_size_test = int(sample_size_test/num_batches)
    remainder1 = sample_size_test%batch_size_test
    remainder2 = sample_size_test%num_batches
    remainder = np.max([remainder1, remainder2])
    num_divisible_test = (sample_size_test - remainder)
    batch_indices_test = np.arange(sample_size_test)
    np.random.shuffle(batch_indices_test)
    batch_indices_test = torch.tensor(batch_indices_test).long()
    batch_indices_test = batch_indices_test[:num_divisible_test]

    X_train_batches = X_train[batch_indices_train.reshape(num_batches, -1)]
    X_test_batches = X_test[batch_indices_test.reshape(num_batches, -1)]
    Y_train_batches = Y_train[batch_indices_train.reshape(num_batches, -1)]
    Y_test_batches = Y_test[batch_indices_test.reshape(num_batches, -1)]

    data = [X_train_batches, X_test_batches, Y_train_batches, Y_test_batches]
    return(data, num_batches)

Y_df = pd.read_csv("y.txt", delimiter = "\t", header = 0)
ICD_cols = Y_df.columns.to_numpy()
ICD_cols = ICD_cols[ICD_cols != "eid"]
Y = Y_df[ICD_cols].to_numpy(dtype = int)
out_dim = len(Y[0])

alpha = 1E-3
final_model_r = []
final_model_p = []
name = "final_model_AE_" + str(dop) + "dop_" + str(rs) + "rs_" + str(dim1) + "d1_" + str(dim2) + "d2.txt"
model_path = "AE_final_model_network/" + name[:-4] + ".model"

print("before training the final model (GB):" + str(process.memory_info().rss/1E9))
network = simpleAE(out_dim, [dim1, dim2, 15], dop) 
optimizer = torch.optim.Adam(network.parameters(), lr = 0.0004, weight_decay = 0)
def factor(k): return(1/(k + 1)) 
scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, factor)
batch_size = 500

loss_func = torch.nn.BCELoss()

loss_vals_final_model = []
loss_r_final_model = []
model_path = "AE_final_model_network/" + name[:-4] + ".model"
if os.path.exists(model_path):
    network.load_state_dict(torch.load(model_path))
for k in range(500):
    time_sum = 0
    batched_data, num_batches = get_batches(torch.tensor(Y).float(), 
                                            torch.tensor(Y).float(), 
                                            torch.tensor(Y), 
                                            torch.tensor(Y), 
                                            batch_size)

    #print("pre gradient descent (GB):" + str(process.memory_info().rss/1E9))
    for i in tqdm(range(num_batches)):

        optimizer.zero_grad()

        Xb_final_model = batched_data[0][i]
        Yb_final_model = batched_data[2][i]

        void, Yb_final_model_pred = network(Xb_final_model)
        final_model_loss = loss_func(Yb_final_model_pred, Yb_final_model.reshape(-1, out_dim).float())

        if i%20 == 0:
            loss_vals_final_model.append(final_model_loss.item())                
            r, p = pearsonr(Yb_final_model_pred.reshape(-1).detach().numpy(), Yb_final_model.reshape(-1).detach().numpy())
            loss_r_final_model.append(r)

        if i%20 == 0:
            network.eval()
            void, Yb_final_model_pred1 = network(Xb_final_model)
            final_model_loss_pre_training = loss_func(Yb_final_model_pred1, Yb_final_model.reshape(-1, out_dim).float()).item()
            network.train()

        final_model_loss.backward()
        optimizer.step()

        if i%20 == 0:
            network.eval()
            void, Yb_final_model_pred2 = network(Xb_final_model)
            final_model_loss_post_training = loss_func(Yb_final_model_pred2, Yb_final_model.reshape(-1, out_dim).float()).item()
            if final_model_loss_pre_training <= final_model_loss_post_training:
                print("The learning rate decreased.")
                scheduler.step()
            network.train()

network.eval()
X1 = torch.tensor(Y).float()
Z, Y_rec = network(X1)
Y_rec = Y_rec.detach().numpy()
Y_pred = Y_rec.reshape(-1)
Y1 = torch.tensor(Y).float()
r, p = pearsonr(Y_pred, Y1.reshape(-1).detach().numpy())
print("training r: " + str(r))
final_model_r.append(r)
final_model_p.append(p)
torch.save(network.state_dict(), model_path)

Y_c = Z.detach().numpy()
any_HF_index = np.where(Y_df.columns[1:].to_numpy()=="any_HF")[0][0]
HF_residuals = Y[:, any_HF_index] - Y_rec[:, any_HF_index]
Y_c_res = np.concatenate([Y_c, HF_residuals.reshape(-1, 1)], axis = 1)
PCs2 = np.concatenate((PCs, np.ones((len(PCs), 1))), axis = 1)
PCsT_PCs = np.matmul(PCs2.T, PCs2)
correction = np.eye(len(PCsT_PCs))*np.min(PCsT_PCs)*1E-6
PCsT_PCs_inv = np.linalg.inv(PCsT_PCs + correction)
coef = np.matmul(PCsT_PCs_inv, PCs2.T)
weights = np.matmul(coef, Y_c_res)
Y_est = np.matmul(PCs2, weights)
phenotypes = Y_c_res - Y_est
pheno_path = "phenotypes_NN_R" + str(r) + ".txt"
phenotypes_df = pd.DataFrame(phenotypes)
phenotypes_df["eid"] = Y_df["eid"]
cnames = phenotypes_df.columns.to_numpy() 
phenotypes_df = phenotypes_df[["eid"] + cnames[cnames != "eid"].tolist()]
phenotypes_df.to_csv(pheno_path , sep = "\t", header = True, index = False)

final_model_errs = pd.DataFrame(np.array([loss_vals_final_model, loss_r_final_model]).T)
final_model_errs.columns = ["final_model_loss", "final_model_r"]
final_model_errs.to_csv("AE_final_model_error_df/" + name, sep = "\t", header = True, index = False)

final_model_info = pd.DataFrame(np.array([final_model_r, final_model_p]).T)
final_model_info.columns = ["final_model r", "final_model p"]
final_model_info.to_csv(name, sep = "\t", header = True, index = False)