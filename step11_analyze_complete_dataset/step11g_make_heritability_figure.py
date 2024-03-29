import numpy as np
import pandas as pd
import os
from scipy.stats import wilcoxon
from matplotlib import pyplot as plt
import pdb

def bootstrap(X, B):
    X_b = np.random.choice(X, (B, len(X)))
    mu_b = np.mean(X_b, axis = 1)
    return(np.percentile(mu_b, [2.5, 97.5]))

NN_paths = os.listdir("../step10_get_significant_SNPs_NN/CV_params")
NN_paths = ["../step10_get_significant_SNPs_NN/CV_params" + "/" + path for path in NN_paths if not "nested" in path]
PCA_paths = os.listdir("../step10_get_significant_SNPs_PCA/CV_params")
PCA_paths = ["../step10_get_significant_SNPs_PCA/CV_params" + "/" + path for path in PCA_paths if not "nested" in path]
logistic_PCA_paths = os.listdir("../step10_get_significant_SNPs_logistic_PCA/CV_params")
logistic_PCA_paths = ["../step10_get_significant_SNPs_logistic_PCA/CV_params" + "/" + path for path in logistic_PCA_paths if not "nested" in path]

r_vals_CV = []
best_params = []
for paths in [NN_paths, PCA_paths, logistic_PCA_paths]:
    files = [pd.read_csv(path, delimiter = "\t") for path in paths]
    GBC_r2_all = np.array([df["R squared (GBC)"].to_numpy() for df in files])
    param_means = np.mean(GBC_r2_all, axis = 0)
    best_row = np.argmax(param_means)
    GBC_r2 = np.array([df["R squared (GBC)"].to_numpy()[best_row] for df in files])
    r_vals_CV.append(GBC_r2)
    best_params.append(files[0][["learning rate", "max depth"]].to_numpy()[best_row])

best_params = pd.DataFrame(np.array(best_params))
best_params["algorithm"] = ["autoencoder", "PCA", "logistic PCA"] 
best_params.columns = ["learning rate", "max depth", "algorithm"]
best_params.to_csv("best_GBC_params.txt", sep = "\t", header = True, index = False)

NN_GBC_CV = r_vals_CV[0]
PCA_GBC_CV = r_vals_CV[1]
logistic_PCA_GBC_CV = r_vals_CV[2]

NN_paths = os.listdir("../step10_get_significant_SNPs_NN/CV_params")
NN_paths = ["../step10_get_significant_SNPs_NN/CV_params" + "/" + path for path in NN_paths if "nested" in path]
PCA_paths = os.listdir("../step10_get_significant_SNPs_PCA/CV_params")
PCA_paths = ["../step10_get_significant_SNPs_PCA/CV_params" + "/" + path for path in PCA_paths if "nested" in path]
logistic_PCA_paths = os.listdir("../step10_get_significant_SNPs_logistic_PCA/CV_params")
logistic_PCA_paths = ["../step10_get_significant_SNPs_logistic_PCA/CV_params" + "/" + path for path in logistic_PCA_paths if "nested" in path]

GBC_r_vals_CV_nested = []
LR_r_vals = []
hyperparameters = []
for paths in [NN_paths, PCA_paths, logistic_PCA_paths]:
    files = [pd.read_csv(path, delimiter = "\t") for path in paths]
    GBC_r_vals_CV_nested.append(np.array([df["R squared (GBC nested CV)"].to_numpy()[0] for df in files]))
    LR_r_vals.append(np.array([df["R squared (LR nested CV)"].to_numpy()[0] for df in files]))
    hyperparameters.append(np.array([df[['learning rate', 'max depth',]].to_numpy()[0] for df in files]))

df_lr, df_md = np.array(hyperparameters).transpose([2, 1, 0])
lr, md = np.arange(0.01, 0.1, 0.01), np.arange(1, 8)
lr_counts = pd.DataFrame(np.sum(df_lr == lr.reshape(-1, 1, 1), axis = 1))
lr_counts.columns = ["autoencoder", "PCA", "logistic PCA"]
lr_counts["learning rate"] = np.round(lr, 2)
lr_counts.to_csv("tableS5b.txt", sep = "\t", header = True, index = False)
lr_counts = lr_counts[["learning rate", "PCA", "logistic PCA", "autoencoder"]]
md_counts = pd.DataFrame(np.sum(df_md == md.reshape(-1, 1, 1), axis = 1))
md_counts.columns = ["autoencoder", "PCA", "logistic PCA"]
md_counts["max depth"] = md
md_counts = md_counts[["max depth", "PCA", "logistic PCA", "autoencoder"]]
md_counts.to_csv("tableS5a.txt", sep = "\t", header = True, index = False)

NN_GBC = GBC_r_vals_CV_nested[0]
NN_LR = LR_r_vals[0]
PCA_GBC = GBC_r_vals_CV_nested[1]
PCA_LR = LR_r_vals[1]
logistic_PCA_GBC = GBC_r_vals_CV_nested[2]
logistic_PCA_LR = LR_r_vals[2]

p_diff_NN = wilcoxon(NN_GBC - NN_LR , alternative = 'greater')[1]
p_diff_PCA = wilcoxon(PCA_GBC - PCA_LR , alternative = 'greater')[1]
p_diff_logistic_PCA = wilcoxon(logistic_PCA_GBC - logistic_PCA_LR, alternative = 'greater')[1]
pvals = np.round([p_diff_NN, p_diff_PCA, p_diff_logistic_PCA], 8).astype(str)

data = [PCA_LR, PCA_GBC, logistic_PCA_LR, logistic_PCA_GBC, NN_LR, NN_GBC]
means = [np.mean(i) for i in data]
inds = np.arange(6) + 1
names = np.array(["PCA (LR)", "PCA (GBC)", "logistic PCA (LR)", "logistic PCA (GBC)", "autoencoder (LR)", "autoencoder (GBC)"])
errs = np.array([bootstrap(i, 100000) for i in data])
max_vals = errs[:, 1][[1, 3, 5]]
errs = np.abs(errs.T - means)

def plot_p_vals(means, inds, pvals, max_vals):
    for k, i in enumerate(range(int(len(means)/2))): 
        inds_sub = inds[2*i:2*(i + 1)]
        max_mean = np.array(2*[np.max(means[2*i:2*(i + 1)])])
        plt.plot(inds_sub, 2*[max_vals[k] + 0.003], "-k", linewidth = 4)
        pval = str(pvals[i])
        if len(pval) == 8: pval = pval[0] + pval[1] + pval[2] + pval[4] + pval[5] + pval[6] + pval[7]
        plt.text(x = inds_sub[0] + 0.08, y = max_vals[k] + 0.006, s = "p=" + pval, fontsize = 24)

fig = plt.figure(figsize=(20, 16))
plt.plot(inds, means, 'ko', markersize = 20)
plt.errorbar(inds, means, yerr = errs, fmt='ko', capsize = 40, elinewidth = 4, capthick = 4)
plt.xticks(inds, labels = names, rotation = 45, ha="right", fontsize=24)
plot_p_vals(means, inds, pvals, max_vals)
title = "within-model mean R squared CIs\n "
title += "with between-model wilcoxon p values"
plt.title(title, size = 40, pad = 20)
plt.tick_params(labelsize = 30)
fig.subplots_adjust(left = 0.15, bottom = 0.3, top = 0.85)
plt.ylabel("R squared between predicted\nand actual heart failure status", fontsize = 32)
plt.ylim([0.04, 0.14])
plt.xlim([0, 7])
plt.savefig("r_squared.png")
plt.clf()
