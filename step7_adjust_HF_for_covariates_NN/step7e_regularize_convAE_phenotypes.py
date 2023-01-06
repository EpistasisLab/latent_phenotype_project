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
    def __init__(self, input_shape, output_shape, latent_dims, dop):
        super(simpleAE, self).__init__()

        dims_in = [output_shape] + latent_dims
        dims_out = np.flip(dims_in)
        self.in_layers = nn.Sequential()
        self.out_layers = nn.Sequential()
        self.in_layers.add_module("first linear", nn.Linear(input_shape, output_shape))
        
        for i in range(len(dims_in) - 1):
             self.in_layers.add_module("in_relu" + str(i), torch.nn.LeakyReLU(negative_slope = 0.05))
             if i > 0: self.in_layers.add_module("in_dropout" + str(i), nn.Dropout(p = dop))
             self.in_layers.add_module("in_linear" + str(i + 1), nn.Linear(dims_in[i], dims_in[i + 1]))
        for i in range(len(dims_in) - 1):
             layer_i = nn.Linear(dims_out[i], dims_out[i + 1])
             layer_i.weight.data = self.in_layers[-3*(i + 1) + 2].weight.data.transpose(0,1)
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

alpha = 1E-3
train_r = []
train_p = []
test_r = []
test_p = []
training_inds = pd.read_csv("training_inds.txt", delimiter = "\t").to_numpy().T
testing_inds = pd.read_csv("testing_inds.txt", delimiter = "\t").to_numpy().T
name = "CV_" + str(dop) + "dop_" + str(rs) + "rs_"  + str(ws) + "ws_" + str(dim1) + "d1_" + str(dim2) + "d2.txt"
model_path = "dropout_CV_fold_networks/" + name[:-4] + ".model"
folds = [1, 2, 3, 4, 5]
for train_inds, test_inds, CV in zip(training_inds, testing_inds, folds):

    ICD_seqs_train = [ICD_seqs[i] for i in train_inds]
    ICD_seqs_test = [ICD_seqs[i] for i in test_inds]
    path = "word2vec_window" + str(ws) + "_seed" + str(rs) + "_CV" + str(CV) + ".model"   
    print(path)
    vec_size = 50
    model =  Word2Vec.load(path)    

    if CV == 1:
        embeddings = model.wv[ICD_codes.tolist()]
        U = UMAP(n_neighbors = 5, n_components = 2, random_state = 0)
        Z = U.fit_transform(embeddings)
        I0 = np.array([seq[0:2] == "I0" for seq in ICD_codes])
        I1 = np.array([seq[0:2] == "I1" for seq in ICD_codes])

        I25 = np.array([seq[0:3] == "I20" for seq in ICD_codes])
        I25 += np.array([seq[0:3] == "I21" for seq in ICD_codes])
        I25 += np.array([seq[0:3] == "I22" for seq in ICD_codes])
        I25 += np.array([seq[0:3] == "I23" for seq in ICD_codes])
        I25 += np.array([seq[0:3] == "I24" for seq in ICD_codes])
        I25 += np.array([seq[0:3] == "I25" for seq in ICD_codes])

        I28 = np.array([seq[0:3] == "I26" for seq in ICD_codes])
        I28 += np.array([seq[0:3] == "I27" for seq in ICD_codes])
        I28 += np.array([seq[0:3] == "I28" for seq in ICD_codes])

        I33 = np.array([seq[0:3] == "I30" for seq in ICD_codes])
        I33 += np.array([seq[0:3] == "I31" for seq in ICD_codes])
        I33 += np.array([seq[0:3] == "I32" for seq in ICD_codes])
        I33 += np.array([seq[0:3] == "I33" for seq in ICD_codes])

        I37 = np.array([seq[0:3] == "I34" for seq in ICD_codes])
        I37 += np.array([seq[0:3] == "I35" for seq in ICD_codes])
        I37 += np.array([seq[0:3] == "I36" for seq in ICD_codes])
        I37 += np.array([seq[0:3] == "I37" for seq in ICD_codes])

        I41 = np.array([seq[0:3] == "I38" for seq in ICD_codes])
        I41 += np.array([seq[0:3] == "I39" for seq in ICD_codes])
        I41 += np.array([seq[0:3] == "I40" for seq in ICD_codes])
        I41 += np.array([seq[0:3] == "I41" for seq in ICD_codes])

        I43 = np.array([seq[0:3] == "I42" for seq in ICD_codes])
        I43 += np.array([seq[0:3] == "I43" for seq in ICD_codes])

        I45 = np.array([seq[0:3] == "I44" for seq in ICD_codes])
        I45 += np.array([seq[0:3] == "I45" for seq in ICD_codes])

        I49 = np.array([seq[0:3] == "I46" for seq in ICD_codes])
        I49 += np.array([seq[0:3] == "I47" for seq in ICD_codes])
        I49 += np.array([seq[0:3] == "I48" for seq in ICD_codes])
        I49 += np.array([seq[0:3] == "I49" for seq in ICD_codes])

        I5 = np.array([seq[0:2] == "I5" for seq in ICD_codes])
        I6 = np.array([seq[0:2] == "I6" for seq in ICD_codes])

        I72 = np.array([seq[0:3] == "I70" for seq in ICD_codes])
        I72 += np.array([seq[0:3] == "I71" for seq in ICD_codes])
        I72 += np.array([seq[0:3] == "I72" for seq in ICD_codes])

        I79 = np.array([seq[0:3] == "I73" for seq in ICD_codes])
        I79 += np.array([seq[0:3] == "I74" for seq in ICD_codes])
        I79 += np.array([seq[0:3] == "I75" for seq in ICD_codes])
        I79 += np.array([seq[0:3] == "I76" for seq in ICD_codes])
        I79 += np.array([seq[0:3] == "I77" for seq in ICD_codes])
        I79 += np.array([seq[0:3] == "I78" for seq in ICD_codes])
        I79 += np.array([seq[0:3] == "I79" for seq in ICD_codes])

        I8 = np.array([seq[0:2] == "I8" for seq in ICD_codes])
        I9 = np.array([seq[0:2] == "I9" for seq in ICD_codes])
        # JR0 = np.array([seq[0:1] == "J" or seq[0:2] == "R0"for seq in ICD_codes])
        # R1_9 = np.array([seq[0:1] == "R" and seq[0:2] != "R0"for seq in ICD_codes])
        groups = [I0, I1, I25, I28, I33, I37, I41, I43, I45, I49, I5, I6, I72, I79, I8, I9]
        marks = ["go", "ro", "bo", "ko", "co", "yo", "gx", "rx", "bx", "kx", "cx", "yx", "g^", "r^", "b^", "k^", "c^", "y^"]
        labs = ["I0", "I1", "I25", "I28", "I33", "I37"]
        labs += ["I41", "I43", "I45", "I49", "I5", "I6", "I72", "I79"]
        labs += ["I8", "I9"]

    X = np.zeros((len(ICD_seqs), vec_size), dtype = np.float32)
    for i, seq in enumerate(ICD_seqs): X[i, :] = np.sum(model.wv[seq], axis = 0)
    print("post ICD embedding RAM (GB):" + str(process.memory_info().rss/1E9))
    N, P = len(X), len(X[0])

    #---------------------------------------------------------------------------------------
    # simulated data test with known R^2 = 0.5
    # pdb.set_trace()
    # N_test = (len(train_inds) + len(test_inds))
    # X = np.random.normal(0, 1, (N_test, 50))
    # W1 = np.random.normal(0, 1, (50, 15))
    # Z = np.matmul(X, W1) 
    # W2 = np.random.normal(0, 1, (15, 311))
    # Y = np.matmul(Z, W2) 
    # Y += np.random.normal(0, np.std(Y), (N_test, 311))
    # X = torch.tensor(X).float()
    # Y = torch.tensor(Y).float()
    #---------------------------------------------------------------------------------------

    X_train, X_test = X[train_inds], X[test_inds]
    Y_train, Y_test = Y[train_inds], Y[test_inds]
    print("post train test splitting (GB):" + str(process.memory_info().rss/1E9))

    network = simpleAE(50, out_dim, [dim1, dim2, 15], dop) 
    optimizer = torch.optim.Adam(network.parameters(), lr = 0.0004, weight_decay = 0)
    def factor(k): return(1/(k + 1)) 
    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, factor)
    batch_size = 500

    loss_func_train = torch.nn.BCELoss()
    loss_func_test = torch.nn.BCELoss()
    #loss_func_train = torch.nn.MSELoss()
    #loss_func_test = torch.nn.MSELoss()

    loss_vals_train = []
    loss_vals_test = []
    loss_r_train = []
    loss_r_test = []

    model_path = "dropout_CV_fold_networks/fold" + str(CV) + "_" + name[:-4] + ".model"
    if os.path.exists(model_path):
        network.load_state_dict(torch.load(model_path))
    else:
        for k in range(200):
            time_sum = 0
            batched_data, num_batches = get_batches(torch.tensor(X_train), 
                                                    torch.tensor(X_test), 
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
                void, Yb_test_pred = network(Xb_test)

                train_loss = loss_func_train(Yb_train_pred, Yb_train.reshape(-1, out_dim).float())
                test_loss = loss_func_test(Yb_test_pred, Yb_test.reshape(-1, out_dim).float())

                loss_vals_train.append(train_loss.item())
                loss_vals_test.append(test_loss.item())
                if i%20 == 0:
                    r, p = pearsonr(Yb_train_pred.reshape(-1).detach().numpy(), Yb_train.reshape(-1).detach().numpy())
                    loss_r_train.append(r)
                    r, p = pearsonr(Yb_test_pred.reshape(-1).detach().numpy(), Yb_test.reshape(-1).detach().numpy())
                    loss_r_test.append(r)

                else:
                    loss_r_train.append(loss_r_train[-1])
                    loss_r_test.append(loss_r_test[-1])
                #start = time()
                train_loss.backward()
                #time_sum += time() - start
                optimizer.step()

                if i%20 == 0:
                    void, Yb_train_pred2 = network(Xb_train)
                    void, Yb_test_pred2 = network(Xb_test)
                    train_loss_post_training = loss_func_train(Yb_train_pred2, Yb_train.reshape(-1, out_dim).float())
                    if train_loss.item() <= train_loss_post_training.item():
                        # print("The learning rate decreased.")
                        scheduler.step()

            #print("gradient descent " + str(i) + " (GB):" + str(process.memory_info().rss/1E9))
            #print(time_sum)

    # note: network.eval() turns off the dropout. network.train() would turn it back on.
    network.eval()
    X1 = torch.tensor(X_train[:int(len(X_test)/1)]).float()
    void, Y_train_pred = network(X1)
    del void
    del X1
    Y_train_pred = Y_train_pred.reshape(-1).detach().numpy()
    Y1 = torch.tensor(Y_train[:int(len(X_test)/1)]).float()
    r, p = pearsonr(Y_train_pred, Y1.reshape(-1).detach().numpy())
    print("training r: " + str(r))
    train_r.append(r)
    train_p.append(p)
    del Y1

    X2 = torch.tensor(X_test[:int(len(X_test)/1)]).float()  
    void, Y_test_pred = network(X2)
    del void
    del X2
    Y_test_pred = Y_test_pred.reshape(-1).detach().numpy()
    Y2 = torch.tensor(Y_test[:int(len(X_test)/1)]).float()
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
    fold_errs.to_csv("dropout_fold_error_dfs/fold" + str(CV) + "_" + name, sep = "\t", header = True, index = False)

CV_info = pd.DataFrame(np.array([train_r, test_r, train_p, test_p]).T)
CV_info.columns = ["training r", "testing r", "training p", "testing p"]
CV_info.to_csv(name, sep = "\t", header = True, index = False)

embeddings2, void = network(torch.tensor(embeddings).float())
embeddings2 = embeddings2.detach().numpy()
U = UMAP(n_neighbors = 5, n_components = 2, random_state = 0)
Z = U.fit_transform(embeddings2)
for g, m, lab in zip(groups, marks, labs): plt.plot(Z[g, 0], Z[g, 1], m, label = lab)
plt.legend()
image_name = "aaa" + str(dop) + "dop_" + str(rs) + "rs_"  + str(ws) + "ws_" + str(dim1) + "d1_" + str(dim2) + "d2.png"
plt.savefig(image_name)
plt.clf()
