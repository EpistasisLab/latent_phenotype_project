import numpy as np
import pandas as pd
import os 
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy.stats import beta
from scipy.stats import chi2
from scipy.stats import wilcoxon
from matplotlib import pyplot as plt
import pdb

all_TRACE_GxE_pvals = []

for k in range(5): 

    fnames = np.array(["step0_GxE_p_values/" + i for i in os.listdir("step0_GxE_p_values") if i[0] == str(k)])
    fname_inds = np.array([int(i.split(".")[0].split("/")[1].split("_")[1]) for i in fnames])
    fnames = fnames[np.argsort(fname_inds)]
    files = [pd.read_csv(i, delimiter = "\t", header = None) for i in fnames]
    TRACE_pvals = np.array([i.to_numpy()[0][0] for i in files])
    all_TRACE_GxE_pvals.append(TRACE_pvals[100:200])
    TRACE_GxE_null_pvals = TRACE_pvals[200:]

    N = len(TRACE_GxE_null_pvals)
    medians = beta.ppf(0.5, np.arange(1, N+1), N - np.arange(N))
    medians = -np.log10(np.flip(medians))
    ubs = beta.ppf(0.975, np.arange(1, N+1), N - np.arange(N))
    ubs = -np.log10(np.flip(ubs))
    lbs = beta.ppf(0.025, np.arange(1, N+1), N - np.arange(N))
    lbs = -np.log10(np.flip(lbs))

    pvals = -np.log10(np.flip(np.sort(TRACE_GxE_null_pvals)))
    chi_alt = chi2.isf(10**(-np.median(pvals)), 1, loc=0, scale=1)
    chi_null = chi2.isf(10**(-np.median(medians)), 1, loc=0, scale=1)
    IF = chi_alt/chi_null

    fig, ax = plt.subplots(figsize = (10, 7.5))
    plt.title('TRACE GxE QQ plot, LP'+ str(k) + ' (IF = ' + str(IF)[0:5] + ')', fontsize=20)

    sorted_p_vals = np.sort(pvals)
    pvals_cdf = np.cumsum(np.ones(len(sorted_p_vals))/len(sorted_p_vals))
    ax.plot(sorted_p_vals, pvals_cdf, "b-")
    ax.plot(sorted_p_vals, np.ones(len(pvals_cdf)), "k-")

    ax2 = ax.twinx()
    ax2.plot(medians, pvals, "g-")
    ax2.plot(medians, lbs, "r-")
    ax2.plot(medians, ubs, "r-")
    ax2.plot(medians, medians, "k-")
    ax2.fill_between(medians, lbs, ubs, color = "r", alpha = 0.5)

    ax.set_xlim([0, np.max(medians)])
    ax.set_xlabel("expected p values under null hypothesis (95% CI)", color = 'r', fontsize=20)
    ax.set_ylabel("CDF of actual p-values", color = 'b', fontsize=20)
    ax2.set_ylabel("actual p-values (QQ plot)", color = 'g', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize = 20)
    ax2.tick_params(axis='both', which='major', labelsize = 20)
    plt.tight_layout()
    plt.savefig("step12c_QQ_plot_TRACE_GxE_null_LP" + str(k) + ".png")
    plt.clf()

all_TRACE_GxE_pvals = np.array(all_TRACE_GxE_pvals)
best_TRACE_pvals = -np.log10(np.min(all_TRACE_GxE_pvals, axis = 0))
all_normal_GxE_pvals = pd.read_csv("step12a_normal_GxE_eff_pvals.txt", delimiter = "\t", header = None)
best_normal_pvals = -np.log10(np.min(all_normal_GxE_pvals, axis = 1))
p_diff = wilcoxon(best_TRACE_pvals, best_normal_pvals, alternative='two-sided')[1]
p_diff = str(p_diff)[0:4] + "e" + str(p_diff).split("e")[1]

plt.hist(best_TRACE_pvals, bins=20, label="best TRACE GxE p-values")
plt.hist(best_normal_pvals, bins=20, label="best normal GxE p-values")
plt.title('TRACE vs normal GxE p-values (p-diff = ' + p_diff + ')', fontsize=20)
plt.xlabel("-log10(p-value)", fontsize=20)
plt.ylabel("-bin count", fontsize=20)
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.subplots_adjust(bottom=0.1, top = 0.925)
plt.yticks(fontsize=14)
plt.savefig("step12c_TRACE vs normal GxE p-values.png")
plt.clf()
