import numpy as np
import pandas as pd
import os
import pdb
from scipy.stats import spearmanr
from scipy.stats import wilcoxon

CV_files = [f for f in os.listdir() if "CV" in f and ".txt" in f]
CV_data = [pd.read_csv(f, delimiter = "\t") for f in CV_files] 

window_sizes = np.array([int(f.split("ws")[0][-1]) for f in CV_files]*5)
AE_shapes = np.array([int(f.split("d1")[0][-4:]) for f in CV_files]*5)
r_test = np.array([df["testing r"].to_numpy() for df in CV_data]).T.reshape(-1)
spearmanr(AE_shapes, r_test)
spearmanr(window_sizes, r_test)

wilcoxon(r_test[np.logical_or(window_sizes == 3, window_sizes == 4)] - r_test[window_sizes <= 2], alternative = "greater")
pdb.set_trace()
print(1)