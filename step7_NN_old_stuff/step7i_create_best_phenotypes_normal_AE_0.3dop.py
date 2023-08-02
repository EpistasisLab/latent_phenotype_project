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
parser.add_argument('--ws', nargs = 1, type = int, action = "store", dest = "ws")
parser.add_argument('--rs', nargs = 1, type = int, action = "store", dest = "rs")
parser.add_argument('--dop', nargs = 1, type = float, action = "store", dest = "dop")

args = parser.parse_args()
ws = args.ws[0]
rs = args.rs[0]
dim1 = args.dim1[0]
dim2 = int(dim1/2)
dop = args.dop[0]

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

prefix = "../step1_get_phenotypes_simple/"
death_ICD_data_all = pd.read_csv(prefix + "death_ICD_data.txt", delimiter = "\t", low_memory = False)
filtered_eids = pd.read_csv("y.txt", delimiter = "\t", usecols = ["eid"])
death_ICD_data_all = death_ICD_data_all.merge(filtered_eids, on = "eid", how = "inner")
del filtered_eids
data_cols = death_ICD_data_all.columns.to_numpy()[death_ICD_data_all.columns.to_numpy() != "eid"]
death_ICD_data = death_ICD_data_all[data_cols]
del death_ICD_data_all

ICD_seqs = [row[pd.isnull(row) == False].tolist() for row in death_ICD_data.to_numpy()]
ICD_seqs = [row if len(row) > 0 else [""] for row in ICD_seqs]
ICD_codes = death_ICD_data.to_numpy().reshape(-1)
del death_ICD_data
print("post data import usage RAM (GB):" + str(process.memory_info().rss/1E9))
ICD_codes = ICD_codes[pd.isnull(ICD_codes) == False]
ICD_codes = np.unique(ICD_codes)

ICD_codes2 = np.array(ICD_codes).reshape(-1, 1)
Y = np.zeros((len(ICD_seqs), len(ICD_codes)), dtype = bool)
for i, seq in enumerate(ICD_seqs): Y[i, :] = np.any(np.array(seq).reshape(-1, 1) == np.array(ICD_codes), axis = 0)
print("post one hot encoding RAM (GB):" + str(process.memory_info().rss/1E9))
out_dim = len(Y[0])

name = "final_model_normalAE_" + str(dop) + "dop_" + str(rs) + "rs_"  + str(ws) + "ws_" + str(dim1) + "d1_" + str(dim2) + "d2.txt"
model_path = "normalAE_pdb_fold_networks/" + name[:-4] + ".model"
    
print("post train test splitting (GB):" + str(process.memory_info().rss/1E9))

network = simpleAE(out_dim, [dim1, dim2, 15], dop) 
optimizer = torch.optim.Adam(network.parameters(), lr = 0.0004, weight_decay = 0)
def factor(k): return(1/(k + 1)) 
scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, factor)
batch_size = 500

loss_func_train = torch.nn.BCELoss()
loss_func_test = torch.nn.BCELoss()

loss_vals_train = []
loss_r_train = []

model_path = "normalAE_final_model_fold_networks/" + name[:-4] + ".model"
if os.path.exists(model_path):
    network.load_state_dict(torch.load(model_path))
for k in range(500):
    time_sum = 0
    # Training the final model does not carve out test folds
    batched_data, num_batches = get_batches(torch.tensor(Y).float(), 
                                            torch.tensor(Y).float(), 
                                            torch.tensor(Y), 
                                            torch.tensor(Y), 
                                            batch_size)

    #print("pre gradient descent (GB):" + str(process.memory_info().rss/1E9))
    for i in tqdm(range(num_batches)):

        optimizer.zero_grad()

        Xb_train = batched_data[0][i]
        # Xb_test = batched_data[1][i]
        Yb_train = batched_data[2][i]
        # Yb_test = batched_data[3][i]

        void, Yb_train_pred = network(Xb_train)
        train_loss = loss_func_train(Yb_train_pred, Yb_train.reshape(-1, out_dim).float())

        if i%20 == 0:
            loss_vals_train.append(train_loss.item())                    
            r, p = pearsonr(Yb_train_pred.reshape(-1).detach().numpy(), Yb_train.reshape(-1).detach().numpy())
            loss_r_train.append(r)

        if i%20 == 0:
            network.eval()
            void, Yb_train_pred1 = network(Xb_train)
            train_loss_pre_training = loss_func_train(Yb_train_pred1, Yb_train.reshape(-1, out_dim).float()).item()
            network.train()

        #start = time()
        train_loss.backward()
        #time_sum += time() - start
        optimizer.step()

        if i%20 == 0:
            network.eval()
            void, Yb_train_pred2 = network(Xb_train)
            train_loss_post_training = loss_func_train(Yb_train_pred2, Yb_train.reshape(-1, out_dim).float()).item()
            if train_loss_pre_training <= train_loss_post_training:
                print("The learning rate decreased.")
                scheduler.step()
            network.train()

    #print("gradient descent " + str(i) + " (GB):" + str(process.memory_info().rss/1E9))
    #print(time_sum)

# note: network.eval() turns off the dropout. network.train() would turn it back on.
network.eval()
X1 = torch.tensor(Y).float()
latent_phenotypes, Y_pred = network(X1)
del X1
Y_pred = Y_pred.reshape(-1).detach().numpy()
Y1 = torch.tensor(Y).float()
model_r, model_p = pearsonr(Y_pred, Y1.reshape(-1).detach().numpy())
print("model r: " + str(model_r))
del Y1

if not os.path.exists(model_path):
    torch.save(network.state_dict(), model_path)

fold_errs = pd.DataFrame(np.array([loss_vals_train, loss_r_train]).T)
fold_errs.columns = ["model loss", "model r"]
fold_errs.to_csv("normalAE_final_model_error_df/" + name, sep = "\t", header = True, index = False)

latent_phenotypes = pd.DataFrame(latent_phenotypes.detach().numpy())
fname = "step7i_" + str(model_r) + "r_output_latent_phenotypes.txt"
latent_phenotypes.to_csv(fname, sep = "\t", header = True, index = False)