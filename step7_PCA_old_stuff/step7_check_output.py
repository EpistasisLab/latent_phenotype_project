import numpy as np
import pandas as pd
import pdb
import statsmodels.api as sm
from bed_reader import open_bed
from scipy.stats import pearsonr
from tqdm import tqdm

def get_effect(x2, y2, E2):
    X2 = np.array([x2, E2, (x2*E2), np.ones(len(E2))]).T
    nans_X2 = np.any(np.isnan(X2), axis = 1)
    nans_y2 = np.isnan(y2)
    vals = np.logical_or(nans_y2, nans_X2) == False
    X, y = X2[vals], y2[vals]
    basic_model_p = sm.OLS(y, X[:, [0, 3]]).fit().pvalues[0]
    model = sm.OLS(y, X)
    model_sub = sm.OLS(y, X[:, [0,1,3]])
    model_results = model.fit()
    model_sub_results = model_sub.fit()
    pval = model_results.compare_lr_test(model_sub_results)[1]
    rsq = model_results.rsquared
    ydev = np.abs(y - np.median(y))
    p_vQTL = pearsonr(ydev, x2[vals])[1]
    return([basic_model_p, rsq, pval, p_vQTL])


is_male = pd.read_csv("X.txt", delimiter = "\t", usecols = ["eid", "22001-0.0"])
is_male.columns = ["eid", "is_male"]
env_factors = pd.read_csv("env_factors.txt", delimiter = "\t", usecols = ["eid", "33", "34"])
phenotypes = pd.read_csv("phenotypes.txt", delimiter = "\t", usecols = ["eid", "16"])
condition1 = np.all(is_male["eid"] == env_factors["eid"])
condition2 = np.all(env_factors["eid"] == phenotypes["eid"])
if not condition1 and condition2:
    print("exiting: these documents must use the same eids column")
    exit()

prefix = "../step8_get_imputed_ukb_samples/filtered_output/UKB_samples_chr6"
col_indices = np.cumsum(40*[18000])
fam_file = pd.read_csv(prefix + ".fam", delim_whitespace = True, header = None)
bim_file = pd.read_csv(prefix + ".bim", delim_whitespace = True, header = None)
genotypes = open_bed(prefix  + ".bed", count_A1 = False, num_threads = 1).read(np.s_[np.argsort(fam_file[1]), col_indices])
non_mafs = np.where(np.nanmean(genotypes, axis = 0)/2 >= 0.5)
for i in non_mafs: genotypes[:, i] = 2 - genotypes[:, i] 
mafs = np.nanmean(genotypes, axis = 0)/2
good_cols = np.logical_and(mafs >= 0.3, np.sum(np.isnan(genotypes), axis = 0) < 5000)
# for i in range(len(genotypes[0])): genotypes[np.isnan(genotypes[:, i]), i] = np.nanmedian(genotypes[:, i]) 

G1, G2 = genotypes[:, good_cols][:, 0:2].T
E1, E2 = env_factors.loc[:, ["33", "34"]].dropna().to_numpy().T
P = phenotypes["16"].dropna().to_numpy()
males = is_male.loc[np.isnan(phenotypes["16"]) == False, "is_male"].to_numpy(dtype = bool)
females = males == False
m = get_effect(G1[males], P[males], E1[males])[2]
f = get_effect(G1[females], P[females], E1[females])[2]
test_ratio = np.max([np.log10(m)/np.log10(f), np.log10(f)/np.log10(m)])

ratios = []
for i in tqdm(range(1000)):
    np.random.shuffle(males)
    females = males == False
    m = get_effect(G1[males], P[males], E1[males])[2]
    f = get_effect(G1[females], P[females], E1[females])[2]
    ratios.append(np.max([np.log10(m)/np.log10(f), np.log10(f)/np.log10(m)]))

pdb.set_trace()
print(1)