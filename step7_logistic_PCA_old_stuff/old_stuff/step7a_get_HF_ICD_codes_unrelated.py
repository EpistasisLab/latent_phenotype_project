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
