import numpy as np
import pandas as pd
import os
import torch
import torch.nn as nn
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from sklearn.model_selection import KFold
from copy import deepcopy as COPY
from tqdm import tqdm
from matplotlib import pyplot as plt
from gensim.test.utils import common_texts
from gensim.models import Word2Vec
from gensim.models.callbacks import CallbackAny2Vec
from sklearn.model_selection import KFold
from sklearn.cluster import KMeans
from umap import UMAP
import pdb
import sys

torch.manual_seed(0)
np.random.seed(0)

class callback(CallbackAny2Vec):
    """
    Callback to print loss after each epoch
    """
    def __init__(self):
        self.epoch = 0

    def on_epoch_end(self, model):
        loss = model.get_latest_training_loss()
        if self.epoch == 0:
            print('Loss after epoch {}: {}'.format(self.epoch, loss))
        else:
            print('Loss after epoch {}: {}'.format(self.epoch, loss- self.loss_previous_step))
        self.epoch += 1
        self.loss_previous_step = loss

class convAE(torch.nn.Module):

    # latent dims last dimension is the output shape
    def __init__(self, input_shape, output_shape, seq_len, latent_dims, num_filters):
        super(convAE, self).__init__()

        k = 5
        self.conv = nn.Conv2d(in_channels = 1, kernel_size = k, out_channels = num_filters)
        self.pool = torch.nn.MaxPool2d(kernel_size = k, padding = int(k/2))

        self.N, self.P = input_shape
        self.N_out = output_shape
        D1 = int(((self.N - (k - 1) + 2*int(k/2))/k))
        D2 = int(((self.P - (k - 1) + 2*int(k/2))/k))
        L = int(D1*D2*num_filters)

        dims_in = [L] + latent_dims
        dims_out = np.flip([output_shape*seq_len] + latent_dims)
        
        self.in_layers = nn.Sequential()
        self.out_layers = nn.Sequential()
        for i in range(len(dims_in) - 1):
             self.in_layers.add_module("in_relu" + str(i), torch.nn.LeakyReLU(negative_slope = 0.05))
             self.in_layers.add_module("in_linear" + str(i + 1), nn.Linear(dims_in[i], dims_in[i + 1]))
        for i in range(len(dims_out) - 1):
             self.out_layers.add_module("out_relu" + str(i), torch.nn.LeakyReLU(negative_slope = 0.05))
             self.out_layers.add_module("out_linear" + str(i + 1), nn.Linear(dims_out[i], dims_out[i + 1]))

        '''
             # for some tied weights
             if i < len(dims_out) - 2:
                 layer_i = nn.Linear(dims_out[i], dims_out[i + 1])
                 layer_i.weight.data = self.in_layers[-2*(i + 1) + 1].weight.data.transpose(0,1)
                 self.out_layers.add_module("out_relu" + str(i), torch.nn.LeakyReLU(negative_slope = 0.05))
                 self.out_layers.add_module("out_linear" + str(i + 1), layer_i)
             else:
                 self.out_layers.add_module("out_relu" + str(i), torch.nn.LeakyReLU(negative_slope = 0.05))
                 self.out_layers.add_module("out_linear" + str(i + 1), nn.Linear(dims_out[i], dims_out[i + 1]))
        '''

        self.sig = nn.Sigmoid()

    def forward(self, X):

        Z = self.pool(self.conv(X)).view(len(X), -1)
        Z2 = self.in_layers(Z)
        # Y_pred = self.sig(self.out_layers(Z2).reshape(-1, self.N_out))
        Y_pred = self.out_layers(Z2).reshape(-1, self.N_out)
        return(Z2, Y_pred)

class simpleAE(torch.nn.Module):

    # latent dims last dimension is the output shape
    def __init__(self, input_shape, output_shape, latent_dims):
        super(simpleAE, self).__init__()

        dims_in = [output_shape] + latent_dims
        dims_out = np.flip(dims_in)
        self.in_layers = nn.Sequential()
        self.out_layers = nn.Sequential()
        self.in_layers.add_module("first linear", nn.Linear(input_shape, output_shape))
        
        for i in range(len(dims_in) - 1):
             self.in_layers.add_module("in_relu" + str(i), torch.nn.LeakyReLU(negative_slope = 0.05))
             self.in_layers.add_module("in_linear" + str(i + 1), nn.Linear(dims_in[i], dims_in[i + 1]))
        for i in range(len(dims_in) - 1):
             layer_i = nn.Linear(dims_out[i], dims_out[i + 1])
             layer_i.weight.data = self.in_layers[-2*(i + 1) + 1].weight.data.transpose(0,1)
             self.out_layers.add_module("out_relu" + str(i), torch.nn.LeakyReLU(negative_slope = 0.05))
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

death_ICD_data_all = pd.read_csv("death_ICD_data.txt", delimiter = "\t", low_memory = False)
data_cols = death_ICD_data_all.columns.to_numpy()[death_ICD_data_all.columns.to_numpy() != "eid"]
death_ICD_data = death_ICD_data_all[data_cols]

ICD_seqs = [row[pd.isnull(row) == False].tolist() for row in death_ICD_data.to_numpy()]
ICD_seqs = [row if len(row) > 0 else [""] for row in ICD_seqs]
ICD_codes = death_ICD_data.to_numpy().reshape(-1)
ICD_codes = ICD_codes[pd.isnull(ICD_codes) == False]
ICD_codes = np.unique(ICD_codes)

seq_lengths = [len(row) for row in ICD_seqs]
L_max = np.max(seq_lengths)
ICD_codes2 = np.array(ICD_codes).reshape(-1, 1)
ICD_ind_seqs = [np.where(np.any(ICD_codes2 == np.array(seq), axis = 1))[0] for seq in ICD_seqs]
seq_lengths2 = [len(row) for row in ICD_ind_seqs]
one_hot_seqs = np.zeros((len(seq_lengths), L_max, len(ICD_codes)), dtype = np.bool)
for i in tqdm(range(len(one_hot_seqs))): 
    for k in range(seq_lengths2[i]):  
        m = ICD_ind_seqs[i][k]
        one_hot_seqs[i][k][m] = True

Y = torch.tensor(one_hot_seqs)
out_dim = len(Y[0][0])

alpha = 1E-3
fold = 1
kf = KFold(n_splits = 5, shuffle = True)
kf.get_n_splits(ICD_seqs)
for train_inds, test_inds in kf.split(ICD_seqs):

    # TODO
    # make python hash file
    # save model with filehash

    ICD_seqs_train = [ICD_seqs[i] for i in train_inds]
    ICD_seqs_test = [ICD_seqs[i] for i in test_inds]
    ws, rs = 1, 0
    path = "word2vec_window_" + str(ws) + "_seed_" + str(rs) + ".model"
    if not os.path.exists(path):
        model = Word2Vec(sentences = ICD_seqs_train, epochs = 100, vector_size = 50, window = ws, seed = rs, alpha = alpha, min_alpha = 1E-6, negative = 20, sg = 1, min_count = 1, workers = 10, compute_loss = True, callbacks=[callback()])
        model.save(path)
    else:
        model =  Word2Vec.load(path)
    
    embeddings = model.wv[ICD_codes.tolist()]
    U = UMAP(n_neighbors = 5, n_components = 2, random_state = 0)
    Z = U.fit_transform(embeddings)
    I0 = np.array([seq[0:2] == "I0" for seq in ICD_codes])
    I1 = np.array([seq[0:2] == "I1" for seq in ICD_codes])
    I2 = np.array([seq[0:2] == "I2" for seq in ICD_codes])
    I3 = np.array([seq[0:2] == "I3" for seq in ICD_codes])
    I4 = np.array([seq[0:2] == "I4" for seq in ICD_codes])
    I5 = np.array([seq[0:2] == "I5" for seq in ICD_codes])
    I6 = np.array([seq[0:2] == "I6" for seq in ICD_codes])
    I7 = np.array([seq[0:2] == "I7" for seq in ICD_codes])
    I8 = np.array([seq[0:2] == "I8" for seq in ICD_codes])
    I9 = np.array([seq[0:2] == "I9" for seq in ICD_codes])
    JR0 = np.array([seq[0:1] == "J" or seq[0:2] == "R0"for seq in ICD_codes])
    R1_9 = np.array([seq[0:1] == "R" and seq[0:2] != "R0"for seq in ICD_codes])
    groups = [I0, I1, I2, I3, I4, I5, I6, I7, I8, I9, JR0, R1_9]
    marks = ["go", "ro", "bo", "ko", "co", "yo", "gx", "rx", "bx", "kx", "cx", "yx",]
    labs = ["I0", "I1", "I2", "I3"]
    labs += ["I4", "I5", "I6", "I7"]
    labs += ["I8", "I9", "JR0", "R1+"]
    for g, m, lab in zip(groups, marks, labs): plt.plot(Z[g, 0], Z[g, 1], m, label = lab)
    plt.legend()
    plt.savefig("aaa.png")
    plt.clf()

    ICD_seq_lengths = [len(seq) for seq in ICD_seqs]
    L_max = np.max(ICD_seq_lengths)
    ICD_seqs2 = [seq + [""]*(L_max - L) for seq, L in zip(ICD_seqs, ICD_seq_lengths)]
    ICD_embeddings = np.array([model.wv[seq] for seq in ICD_seqs2])
    X = torch.tensor(ICD_embeddings).float()
    N, P = len(X[0]), len(X[0][0])
    # X = X.view(len(X), 1, N, P)

    #---------------------------------------------------------------------------------------
    # simulated data test with known R^2 = 0.5
    # pdb.set_trace()
    # del ICD_embeddings
    # del ICD_seqs2
    # N_test = (len(train_inds) + len(test_inds))
    # X = np.random.normal(0, 1, (N_test, 50))
    # W1 = np.random.normal(0, 1, (50, 15))
    # Z = np.matmul(X, W1) 
    # W2 = np.random.normal(0, 1, (15, 311))
    # Y = np.matmul(Z, W2) 
    # Y += np.random.normal(0, np.std(Y), (N_test, 311))
    # X = torch.tensor(X).float()
    # Y = torch.tensor(Y).float()
    #
    # # only needed if trying convAE, which doesn't really work. 
    # X = X.reshape(-1, 1, 50, 50)
    # Y = Y.reshape(-1, 50, 311)
    # inds = np.arange(N_test/50)
    # np.random.shuffle(inds)
    # train_inds = inds[:80000]
    # test_inds = inds[:20000]
    #---------------------------------------------------------------------------------------

    X_train, X_test = X[train_inds], X[test_inds]
    Y_train, Y_test = Y[train_inds], Y[test_inds]

    loss_func_train = torch.nn.BCELoss()
    loss_func_test = torch.nn.BCELoss()
    #loss_func_train = torch.nn.MSELoss()
    #loss_func_test = torch.nn.MSELoss()

    # network = convAE([N, P], out_dim, 50, [100, 15], 50)
    network = simpleAE(50, out_dim, [250, 125, 15]) 
    optimizer = torch.optim.Adam(network.parameters(), lr = 0.001, weight_decay = 0) 
    epochs = 10
    batch_size = 500

    #---------------------------------------------------------------------------------------
    # overfitting test: permuting X_train should make testing_error represent noise
    # old_shape = X_train.shape
    # indices = np.arange(len(X_train.reshape(-1)))
    # np.random.shuffle(indices)
    # indices = torch.tensor(indices).long()
    # X_train = X_train.reshape(-1)[indices].reshape(old_shape)
    #---------------------------------------------------------------------------------------
    
    batched_data, num_batches = get_batches(X_train, X_test, Y_train, Y_test, batch_size)
    loss_vals_train = []
    loss_vals_test = []

    for k in range(10):
        for i in tqdm(range(num_batches)):

            optimizer.zero_grad()

            Xb_train = batched_data[0][i]
            Xb_test = batched_data[1][i]
            Yb_train = batched_data[2][i]
            Yb_test = batched_data[3][i]

            # void, Yb_train_pred = network(Xb_train)
            # void, Yb_test_pred = network(Xb_test)
            void, Yb_train_pred = network(Xb_train.reshape(-1, 50))
            void, Yb_test_pred = network(Xb_test.reshape(-1, 50))

            train_loss = loss_func_train(Yb_train_pred, Yb_train.reshape(-1, out_dim).float())
            test_loss = loss_func_test(Yb_test_pred, Yb_test.reshape(-1, out_dim).float())

            loss_vals_train.append(train_loss.item())
            loss_vals_test.append(test_loss.item())
            train_loss.backward()
            optimizer.step()

    M = 100
    void, Y_train_pred = network(X_train[:len(X_test)].reshape(-1, 50))
    # void, Y_train_pred = network(X_train[:len(X_test)])
    print(pearsonr(Y_train_pred.reshape(-1).detach().numpy(), Y_train[:len(X_test)].reshape(-1).detach().numpy()))

    void, Y_test_pred = network(X_test.reshape(-1, 50))
    # void, Y_test_pred = network(X_test)
    print(pearsonr(Y_test_pred.reshape(-1).detach().numpy(), Y_test.reshape(-1).detach().numpy()))

    '''
    plt.plot(np.arange(len(loss_vals_train[:M])), loss_vals_train[:M], "-", label = "training")
    plt.plot(np.arange(len(loss_vals_test[:M])), loss_vals_test[:M], "-", label = "testing")
    plt.savefig("aaa" + str(fold) + ".png")
    plt.clf()
    '''
   
    embeddings2, void = network(torch.tensor(embeddings).float())
    embeddings2 = embeddings2.detach().numpy()
    U = UMAP(n_neighbors = 5, n_components = 2, random_state = 0)
    Z = U.fit_transform(embeddings2)
    for g, m, lab in zip(groups, marks, labs): plt.plot(Z[g, 0], Z[g, 1], m, label = lab)
    plt.legend()
    plt.savefig("aaa_new.png")
    plt.clf()

    fold += 1
    pdb.set_trace()

print(1)

'''
            weights_train = torch.zeros(Yb_train.shape)
            weights_test = torch.zeros(Yb_test.shape)
            weights_train[Yb_train] = w1
            weights_train[Yb_train == False] = w0
            weights_test[Yb_test] = w1
            weights_test[Yb_test == False] = w0
''' 
