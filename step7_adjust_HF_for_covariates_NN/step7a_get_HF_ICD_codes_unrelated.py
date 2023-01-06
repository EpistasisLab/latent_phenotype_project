from sklearn.model_selection import KFold
import numpy as np
import pandas as pd
import pdb

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

kf = KFold(n_splits = 5, shuffle = True, random_state = 0)
kf.get_n_splits(X)
CV_ind_sets = [i for i in kf.split(X)]
training_inds = pd.DataFrame(np.array([i[0] for i in CV_ind_sets]).T)
testing_inds = pd.DataFrame(np.array([i[1] for i in CV_ind_sets]).T)
training_inds.to_csv("training_inds.txt", sep = "\t", header = True, index = False)
testing_inds.to_csv("testing_inds.txt", sep = "\t", header = True, index = False)