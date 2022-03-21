import numpy as np
import pandas as pd
from scipy.stats import linregress
from matplotlib import pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
from time import time
import numpy as np
import argparse
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('--phe', nargs = 1, type = str, action = "store", dest = "phe")
args = parser.parse_args()
phe = (args.phe)[0]

def jack_knife(X, y, ind_sub, M):
     params, void, void, void = np.linalg.lstsq(X, y, rcond = None)
     slope_main, intercept_main = params

     with Pool(8) as p:
         params = p.map(regress_subset, ind_sub)
     params = np.array(params)
     slopes, intercepts = params[:, 0], params[:, 1]
     slope_modified = M*slope_main - (M - 1)*np.mean(slopes)
     intercept_modified = M*intercept_main - (M - 1)*np.mean(intercepts)
     return(slope_modified, intercept_modified)

chr = ["1", "2", "3", "4", "5", "6", "7", "8", "9"]
chr += ["10", "11", "12", "13", "14", "15", "16"]
chr += ["17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
info = pd.concat([pd.read_csv("LD_scores/chr_" + str(i) + ".txt", delimiter = "\t") for i in chr])
prefix = "QTL_output/QTL_effects_chr"
info2 = pd.concat([pd.read_csv(prefix + i + "_P" + phe + ".txt", delimiter = "\t") for i in chr])

good_vals = (info2["p_null2"] < 0.05).to_numpy()
N, M = np.mean(info.loc[good_vals, "N"]), len(info2)
y_raw = ((info2.loc[good_vals, "p_null2_beta"]**2)*info.loc[good_vals, "N"]).to_numpy()
x_raw = info.loc[good_vals, "ld_score"].to_numpy()
B = 100
sorted_indices = np.argsort(x_raw)
x_sorted, y_sorted = x_raw[sorted_indices], y_raw[sorted_indices]
num_bins = int(len(x_raw)/B)
real_num = len(x_raw)/B
x_medians = [np.mean(x_sorted[i*B:(i+1)*B]) for i in range(num_bins)]
if real_num > num_bins: x_medians.append(np.median(x_sorted[num_bins*B:]))
y_medians = [np.mean(y_sorted[i*B:(i+1)*B]) for i in range(num_bins)]
if real_num > num_bins: y_medians.append(np.median(y_sorted[num_bins*B:]))
plt.plot(x_medians, y_medians, "*")
plt.savefig("bbb.png")
plt.clf()
pdb.set_trace()
 

slope, intercept, r_val, p_val, void = linregress(x, y)
h2 = slope*M/N

ind = np.arange(len(x))
ind_sub = np.arange(len(x)) # np.random.choice(ind, 100, replace = False)
X = np.concatenate([x.reshape(-1, 1), np.ones((len(x), 1))], axis = 1)

def regress_subset(i):
    ind_jack = ind[ind != i]
    params, void, void, void = np.linalg.lstsq(X[ind_jack], y[ind_jack], rcond = None)
    return(params)

slope_adj, intercept_adj = jack_knife(X, y, ind_sub, M)
h2_adj = slope_adj*(M/N)
params = [[slope, intercept, slope_adj, intercept_adj, M, N, h2, h2_adj]]
params += [[slope, intercept, slope_adj, intercept_adj, M, N, h2, h2_adj]]
params = pd.DataFrame(params)
params.columns = ["slope", "intercept", "slope_adj", "intercept_adj", "M", "N", "h2", "h2_adj"]
params.to_csv("LD_scores/LD_reg_pheno" + phe + ".txt", sep = "\t", header = True, index = False)