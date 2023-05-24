import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
import pdb

eid_path = "../step9_regress_phenotypes_against_SNPs_PCA/QTL_phenotypes_eids.txt"
eids = pd.read_csv(eid_path, delimiter = "\t", header = None)
eids.columns = ["eid"]

path = "../step7_adjust_HF_for_covariates_PCA/y.txt"
data = pd.read_csv(path, delimiter = "\t") 
data = data.merge(eids, on = "eid", how = "inner")
y = data['any_HF'].to_numpy(dtype = int)

main_path = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr1.fam"
new_fam_main = pd.read_csv(main_path, delim_whitespace = True, header = None)[1].to_numpy()
is_in_intersect = np.isin(new_fam_main, eids.to_numpy().reshape(-1))
N_total = np.sum(is_in_intersect)
N_cases = np.sum(y)
N_conts = N_total - N_cases

all_inds = np.arange(N_total)
inds_cases = all_inds[y == 1]
inds_conts = all_inds[y == 0] 

kf_cases = KFold(n_splits = 30, shuffle = True)
kf_conts = KFold(n_splits = 30, shuffle = True)
train_inds, test_inds = [], []
for N, kf, inds in zip([N_cases, N_conts], [kf_cases, kf_conts], [inds_cases, inds_conts]):
    kf.get_n_splits(np.ones((N, 2)))
    ind_sets = [ind_set for ind_set in kf.split(np.ones((N, 2)))]
    train_inds.append(pd.DataFrame([inds[i[0]] for i in ind_sets]).T)
    test_inds.append(pd.DataFrame([inds[i[1]] for i in ind_sets]).T)

z1 = pd.concat(train_inds).to_numpy()
z2 = pd.concat(test_inds).to_numpy()
check1 = np.array([len(np.unique(z1[np.isnan(z1[:, i]) == False, i])) + len(np.unique(z2[np.isnan(z2[:, i]) == False, i])) for i in range(30)])
check2 = len(np.unique(z1.reshape(-1)[np.isnan(z1.reshape(-1)) == False]))
if not (np.all(check1 == N_total) and check2 == N_total):
    print("there is a bug in the cross validation index division")
    pdb.set_trace()

pd.concat(train_inds).to_csv("step10f_train_inds.txt", sep = "\t", header = True, index = False)
pd.concat(test_inds).to_csv("step10f_test_inds.txt", sep = "\t", header = True, index = False)




