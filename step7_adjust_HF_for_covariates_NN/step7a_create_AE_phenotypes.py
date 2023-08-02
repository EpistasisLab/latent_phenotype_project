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
train_r = []
train_p = []
test_r = []
test_p = []
training_inds = pd.read_csv("training_inds.txt", delimiter = "\t").to_numpy().T
testing_inds = pd.read_csv("testing_inds.txt", delimiter = "\t").to_numpy().T
name = "CV_AE_" + str(dop) + "dop_" + str(rs) + "rs_" + str(dim1) + "d1_" + str(dim2) + "d2.txt"
model_path = "AE_CV_fold_networks/" + name[:-4] + ".model"
folds = [1, 2, 3, 4, 5]
for train_inds, test_inds, CV in zip(training_inds, testing_inds, folds):

    Y_train, Y_test = Y[train_inds], Y[test_inds]
    print("post train test splitting (GB):" + str(process.memory_info().rss/1E9))
    network = simpleAE(out_dim, [dim1, dim2, 15], dop) 
    optimizer = torch.optim.Adam(network.parameters(), lr = 0.0004, weight_decay = 0)
    def factor(k): return(1/(k + 1)) 
    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, factor)
    batch_size = 500

    loss_func_train = torch.nn.BCELoss()
    loss_func_test = torch.nn.BCELoss()

    loss_vals_train = []
    loss_vals_test = []
    loss_r_train = []
    loss_r_test = []

    model_path = "AE_CV_fold_networks/fold" + str(CV) + "_" + name[:-4] + ".model"
    if os.path.exists(model_path):
        network.load_state_dict(torch.load(model_path))
    for k in range(500):
        time_sum = 0
        batched_data, num_batches = get_batches(torch.tensor(Y_train).float(), 
                                                torch.tensor(Y_test).float(), 
                                                torch.tensor(Y_train), 
                                                torch.tensor(Y_test), 
                                                batch_size)

        #print("pre gradient descent (GB):" + str(process.memory_info().rss/1E9))
        for i in tqdm(range(num_batches)):

            optimizer.zero_grad()

            Xb_train = batched_data[0][i]
            Xb_test = batched_data[1][i]
            Yb_train = batched_data[2][i]
            Yb_test = batched_data[3][i]

            void, Yb_train_pred = network(Xb_train)
            train_loss = loss_func_train(Yb_train_pred, Yb_train.reshape(-1, out_dim).float())

            if i%20 == 0:
                loss_vals_train.append(train_loss.item())

                void, Yb_test_pred = network(Xb_test)
                test_loss = loss_func_test(Yb_test_pred, Yb_test.reshape(-1, out_dim).float())
                loss_vals_test.append(test_loss.item())
                    
                r, p = pearsonr(Yb_train_pred.reshape(-1).detach().numpy(), Yb_train.reshape(-1).detach().numpy())
                loss_r_train.append(r)
                r, p = pearsonr(Yb_test_pred.reshape(-1).detach().numpy(), Yb_test.reshape(-1).detach().numpy())
                loss_r_test.append(r)

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
    X1 = torch.tensor(Y_train[:int(len(Y_test)/1)]).float()
    void, Y_train_pred = network(X1)
    del void
    del X1
    Y_train_pred = Y_train_pred.reshape(-1).detach().numpy()
    Y1 = torch.tensor(Y_train[:int(len(Y_test)/1)]).float()
    r, p = pearsonr(Y_train_pred, Y1.reshape(-1).detach().numpy())
    print("training r: " + str(r))
    train_r.append(r)
    train_p.append(p)
    del Y1

    X2 = torch.tensor(Y_test[:int(len(Y_test)/1)]).float()  
    void, Y_test_pred = network(X2)
    del void
    del X2
    Y_test_pred = Y_test_pred.reshape(-1).detach().numpy()
    Y2 = torch.tensor(Y_test[:int(len(Y_test)/1)]).float()
    r, p = pearsonr(Y_test_pred, Y2.reshape(-1).detach().numpy())
    print("testing r: " + str(r))
    test_r.append(r)
    test_p.append(p)
    del Y2

    '''
    M = 100
    plt.plot(np.arange(len(loss_vals_train[:M])), loss_vals_train[:M], "-", label = "training")
    plt.plot(np.arange(len(loss_vals_test[:M])), loss_vals_test[:M], "-", label = "testing")
    plt.savefig("aaa" + str(CV) + ".png")
    plt.clf()
    '''

    if not os.path.exists(model_path):
        torch.save(network.state_dict(), model_path)

    fold_errs = pd.DataFrame(np.array([loss_vals_train, loss_vals_test, loss_r_train, loss_r_test]).T)
    fold_errs.columns = ["training", "testing", "training r", "testing r"]
    fold_errs.to_csv("AE_fold_error_dfs/fold" + str(CV) + "_" + name, sep = "\t", header = True, index = False)

CV_info = pd.DataFrame(np.array([train_r, test_r, train_p, test_p]).T)
CV_info.columns = ["training r", "testing r", "training p", "testing p"]
CV_info.to_csv(name, sep = "\t", header = True, index = False)