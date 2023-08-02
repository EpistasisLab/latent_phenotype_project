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
parser.add_argument('--ws', nargs = 1, type = int, action = "store", dest = "ws")
parser.add_argument('--rs', nargs = 1, type = int, action = "store", dest = "rs")

args = parser.parse_args()
ws = args.ws[0]
rs = args.rs[0]

# prints ram usage (GB)
process = psutil.Process(os.getpid())
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
fold = 1
train_r = []
train_p = []
test_r = []
test_p = []
training_inds = pd.read_csv("training_inds.txt", delimiter = "\t").to_numpy().T
testing_inds = pd.read_csv("testing_inds.txt", delimiter = "\t").to_numpy().T
folds = [1, 2, 3, 4, 5]
for train_inds, test_inds, CV in zip(training_inds, testing_inds, folds):

    ICD_seqs_train = [ICD_seqs[i] for i in train_inds]
    ICD_seqs_test = [ICD_seqs[i] for i in test_inds]
    path = "word2vec_window" + str(ws) + "_seed" + str(rs) + "_CV" + str(CV) + ".model"
    vec_size = 50
    model = Word2Vec(sentences = ICD_seqs_train, epochs = 1000, vector_size = vec_size, window = ws, seed = rs, alpha = alpha, min_alpha = 1E-6, negative = 20, sg = 1, min_count = 1, workers = 10, compute_loss = True, callbacks=[callback()])
    model.save(path)
    
    if CV == 1:
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
        plt.savefig("word2vec_UMAP_window" + str(ws) + "_seed" + str(rs) + "_CV" + str(CV) + ".png")
        plt.clf()

