import numpy as np
import pandas as pd
import os
import argparse
from copy import deepcopy as COPY
import statsmodels.api as sm
from sklearn.decomposition import PCA
from scipy.linalg import sqrtm
from scipy.stats import pearsonr
from scipy.stats import yeojohnson as yj
from scipy.stats import linregress
from matplotlib import pyplot as plt
from scipy.stats import binom_test
from scipy.stats import chi2
from scipy.stats import rankdata
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import TruncatedSVD as SVD
from tqdm import tqdm
from time import time
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('--k', nargs = 1, type = str, action = "store", dest = "k")
args = parser.parse_args()
k = int((args.k)[0])

def link_pairs(pairs):
#----------------------------------------------------------------------------------
# The goal is to group ICD codes so that all of the covariates in
# each set are at least moderately correlated (|r| >= 0.4) to 
# at least one other covariate within the set, and that all sets
# contain at least two covariates. 
#
# INPUT (pairs): all pairs of covariate indices with |r| >= 0.4
#                all such pairs occur twice as [i, j] and [j, i]
#                all pairs [a, b] are sorted by a and then b
#
# OUTPUT (groups): sets of covariate indices as described above
#
# NOTE: Consider the pair union of index val. Suppose that the pair union 
#       of index (val + x) does not contain val but does contain (val + x + y). 
#       Then suppose the pair union of index (val + x + y) contains val. The
#       concern of forming two seperate groups that should be linked 
#       ([(val + x), (val + x + y)] and [(val + x + y), (val)]) is impossible.
#       If the pair union of index val + x + y contains val, then the
#       pair union of index val must also contain index val + x + y. The pair
#       union of index (val + x) would then be linked to that of val by 
#       (val + x + y), which would form [val, (val + x), (val + x + y)].
#----------------------------------------------------------------------------------
    vals = np.unique(pairs)
    
    groups = []
    group_lengths = np.array([])
    for val in vals: 

        # Only groups with at least one element of the next_union_of_pairs
        # should be unioned with that union of pairs. This is done for each set
        # of pairs with every possible index (val) being the pairs first element
        # in increasing order. This prevents early formulation of new sets (see NOTE)
        next_union_of_pairs = np.unique(pairs[pairs[:, 0] == val, :])
        group_pair_unions = [np.union1d(next_union_of_pairs, g) for g in groups]

        # the union length = pair length + group length if the group
        # shares no elements with the pair. Such unions are omitted. 
        union_lengths = np.array([len(u) for u in group_pair_unions])
        good_union = union_lengths < (group_lengths + len(next_union_of_pairs))

        # The exception is if no groups contain at least one element
        # of the pair. That becomes a new group. 
        if np.all(union_lengths == (group_lengths + len(next_union_of_pairs))):
            groups.append(next_union_of_pairs)
            group_lengths = np.array([len(g) for g in groups])
        else:
            for i in range(len(groups)):
                if good_union[i]:
                    groups[i] = group_pair_unions[i]
            group_lengths = np.array([len(g) for g in groups])
    return(groups)

def get_group_importance(g_ind, g_labels, cors2, rots_original, labels):

    # TODO 1) abandon this function
    #      2) plot hist of top ten most emphasized ICDs
    #      3) stacked bar chart of each most emphasized component in new phenotype
    #         other components not most emphasized can be in "other"

    rots = COPY(rots_original)
    group_importance = []
    for ind in g_ind:
        weights = np.mean(cors2[ind][:, ind], axis = 0)
        group_importance.append(np.sum(weights*rots[ind]))
        rots[ind] *= (1 - weights)

    importance = np.concatenate([rots, group_importance])
    last_labels = np.array(labels.to_list() + g_labels)
    return(importance, last_labels)

def scale_color(vals):
    r = np.max(np.array([-vals + 1, np.ones(len(vals))]), axis = 0)
    r[r > 1] -= 1
    g = np.max(np.array([vals, np.zeros(len(vals))]), axis = 0)
    b = np.ones(len(vals))*(r < 1)
    colors = np.array([r, g, b]).T
    return(colors)

X_df = pd.read_csv("X.txt", delimiter = "\t", header = 0)
X_cols = np.array(X_df.columns)
X_cols[X_cols == "22001-0.0"] = "is_male"
X_df.columns = X_cols
X = X_df[X_cols[X_cols != "eid"]].to_numpy() 
X_cols = X_cols[X_cols != "eid"]
is_PC = np.array([col[0:2] == "PC" for col in X_cols])
PCs = X[:, is_PC]

y_df = pd.read_csv("y.txt", delimiter = "\t", header = 0)
y_cols = y_df.columns
Y = y_df[y_cols[1:]].to_numpy(dtype = float)

# NOTE: all similar files in the logistic_SVD_output folder
#       are made slightly differently, though the results are mostly the same. 
#       these will be easier for others to interpret. 

Y_cors = np.corrcoef(Y.T) 
signs = np.sign(Y_cors)
Y_cors = Y_cors**2 - np.eye(len(Y_cors[0])) 
Y_pairs = np.array(np.where(Y_cors > 0.16)).T
Y_label_sets_ind = link_pairs(Y_pairs)      
Y_label_sets = ["-".join([name.split("s_")[-1] for name in y_cols[1:][ind]]) for ind in Y_label_sets_ind]
Y_cors2 = Y_cors*signs

'''
# raw_X = pd.DataFrame(Y_c)
# raw_Y = pd.DataFrame(Y_rots)
Y_c = pd.read_csv("logistic_SVD_output/raw_X" + str(k) + ".txt", delimiter = "\t", header = None).to_numpy()
Y_rots_inv = pd.read_csv("logistic_SVD_output/raw_Y" + str(k) + ".txt", delimiter = "\t", header = None).to_numpy()
Y_rots = np.linalg.pinv(Y_rots_inv).T
'''

pcaY = (PCA(n_components = k)).fit(Y)
Y_c = pcaY.transform(Y) 
Y_rots = pcaY.components_


'''
# confirm that they are rotations for Y
zz = np.matmul((Y - np.mean(Y, axis = 0)), Y_rots.T).reshape(-1)  
zzz = Y_c.reshape(-1)
pearsonr(zz, zzz)
# output: (0.999999999999994, 0.0)
'''


abs_Y_rots = np.abs(Y_rots)
normalized_Y_rots = abs_Y_rots/np.sum(abs_Y_rots, axis = 1, keepdims = True)
importances = np.sum(normalized_Y_rots, axis = 0)
sorted_indices = np.flip(np.argsort(importances))
sorted_importances = importances[sorted_indices]
sorted_Y_names = y_cols[1:][sorted_indices]
HF_names = np.array(["LV_HF", "CMyo", "any_HF", "Unknown_HF", "cong_HF"])
cutoff = np.max(np.where(np.isin(sorted_Y_names, HF_names))[0])


# real_names correspond to these ICD codes in the listed order
# note: Assuming I7020 is I702
# note: calling "aneurysm of heart" "aortic aneurysm" because that's what it seems to refer to. 
# note: postprocedural disorders = "Postprocedural disorders of circulatory system, not elsewhere classified"
# note: "hypertensive heart and renal disease" is labeled "hypertensive heart disease"
ICD_names = ['counts_I10', 'counts_R074', 'counts_R31', 'counts_R11', 'counts_R69']
ICD_names += ['counts_R55', 'counts_R33', 'counts_I251', 'counts_I48', 'counts_I849']
ICD_names += ['counts_R51', 'counts_I839', 'counts_J90', 'counts_R073', 'counts_R42']
ICD_names += ['counts_I269', 'counts_I959', 'counts_I489', 'counts_I209']
ICD_names += ['counts_I259', 'any_HF', 'counts_R13', 'LV_HF', 'counts_I252']
ICD_names += ['counts_R35', 'counts_R002', 'counts_R060', 'counts_R32', 'counts_R53']
ICD_names += ['counts_I802', 'counts_R05', 'counts_I951', 'counts_R001']
ICD_names += ['counts_I842', 'counts_R15', 'counts_I480', 'counts_I848']
ICD_names += ['counts_I219', 'counts_R91', 'counts_I258', 'counts_I214']
ICD_names += ['counts_I639', 'counts_I200', 'counts_I739', 'counts_R072']
ICD_names += ['counts_R042', 'Unknown_HF', 'counts_R21', 'counts_I841', 'counts_R000']
ICD_names += ['counts_I471', 'cong_HF', 'counts_I350', 'counts_I730', 'CMyo']
ICD_names += ['counts_R040', 'counts_R14', 'counts_R18']
ICD_names += ['counts_R12', 'counts_I64', 'counts_I694'] 
ICD_names += ['counts_I517', 'counts_I211', 'counts_I210']
ICD_names += ['counts_R065', 'counts_I678', 'counts_I440'] 
ICD_names += ['counts_R54', 'counts_I451', 'counts_R011'] 
ICD_names += ['counts_I679', 'counts_I609', 'counts_I831'] 
ICD_names += ['counts_R030', 'counts_I743', 'counts_I340'] 
ICD_names += ['counts_I518', 'counts_I693', 'counts_I249'] 
ICD_names += ['counts_I845', 'counts_I872', 'counts_I7020'] 
ICD_names += ['counts_R17', 'counts_I341', 'counts_I499'] 
ICD_names += ['counts_R02', 'counts_I771', 'counts_I859'] 
ICD_names += ['counts_I472', 'counts_I635', 'counts_I800'] 
ICD_names += ['counts_I120', 'counts_R031', 'counts_I352'] 
ICD_names += ['counts_I447', 'counts_I460', 'counts_R600'] 
ICD_names += ['counts_I801', 'counts_R068', 'counts_I671']
ICD_names += ['counts_I619', 'counts_I208', 'counts_I652']
ICD_names += ['counts_I080', 'counts_I441', 'counts_I890']
ICD_names += ['counts_I702', 'counts_I229', 'counts_I781']
ICD_names += ['counts_I982', 'counts_I248', 'counts_I481']
ICD_names += ['counts_I830', 'counts_I803', 'counts_I493']
ICD_names += ['counts_I469', 'counts_R070', 'counts_I083']
ICD_names += ['counts_I313', 'counts_I719', 'counts_I490']
ICD_names += ['counts_I494', 'counts_R091', 'counts_I456']
ICD_names += ['counts_I850', 'counts_I351', 'counts_J91']
ICD_names += ['counts_R80', 'counts_R161', 'counts_I638']
ICD_names += ['counts_R062', 'counts_I272', 'counts_I442']
ICD_names += ['counts_I371', 'counts_I38', 'counts_I620']
ICD_names += ['counts_I672', 'counts_I779', 'counts_I119']
ICD_names += ['counts_I279', 'counts_I453', 'counts_I212']
ICD_names += ['counts_I712', 'counts_J81', 'counts_I698']
ICD_names += ['counts_I059', 'counts_R093', 'counts_R064']
ICD_names += ['counts_I221', 'counts_I776', 'counts_I495']
ICD_names += ['counts_I7021', 'counts_I714', 'counts_I452']
ICD_names += ['counts_I691', 'counts_I690', 'counts_I213']
ICD_names += ['counts_I809', 'counts_I633', 'counts_I270']
ICD_names += ['counts_I082', 'counts_I864', 'counts_I081']
ICD_names += ['counts_R071', 'counts_I745', 'counts_I519']
ICD_names += ['counts_R098', 'counts_R72', 'counts_I319']
ICD_names += ['counts_I81', 'counts_I482', 'counts_I455']
ICD_names += ['counts_I359', 'counts_R609', 'counts_I479']
ICD_names += ['counts_I634', 'counts_I629', 'counts_R008']
ICD_names += ['counts_I632', 'counts_I443', 'counts_I470']
ICD_names += ['counts_I891', 'counts_I358', 'counts_I483']
ICD_names += ['counts_I822', 'counts_I498', 'counts_I868']
ICD_names += ['counts_R092', 'counts_I446', 'counts_R160']
ICD_names += ['counts_I088', 'counts_I610', 'counts_I618']
ICD_names += ['counts_I878', 'counts_I459', 'counts_I309']
ICD_names += ['counts_I828', 'counts_I952', 'counts_I898']
ICD_names += ['counts_I201', 'counts_I612', 'counts_I7000']
ICD_names += ['counts_I724', 'counts_I832', 'counts_R090']
ICD_names += ['counts_I330', 'counts_I792', 'counts_R061']
ICD_names += ['counts_I983', 'counts_I99', 'counts_I071']
ICD_names += ['counts_I052', 'counts_I958', 'counts_I129']
ICD_names += ['counts_I260', 'counts_I050', 'counts_I780']
ICD_names += ['counts_I511', 'counts_I829', 'counts_I491']
ICD_names += ['counts_R162', 'counts_I254', 'counts_R64']
ICD_names += ['counts_I611', 'counts_I871', 'counts_I748']
ICD_names += ['counts_I630', 'counts_I458', 'counts_I710']
ICD_names += ['counts_I713', 'counts_R58', 'counts_I749']
ICD_names += ['counts_I788', 'counts_I7080', 'counts_I950']
ICD_names += ['counts_I444', 'counts_I220', 'counts_R34']
ICD_names += ['counts_I342', 'counts_I778', 'counts_I253']
ICD_names += ['counts_I978', 'counts_I349', 'counts_I361']
ICD_names += ['counts_I709', 'counts_I808', 'counts_I514']
ICD_names += ['counts_I278', 'counts_I240', 'counts_R048']
ICD_names += ['counts_I723', 'counts_R36', 'counts_I454']
ICD_names += ['counts_I823', 'counts_I744', 'counts_I971']
ICD_names += ['counts_I650', 'counts_R601', 'counts_I728']
ICD_names += ['counts_I692', 'counts_I738', 'counts_I228']
ICD_names += ['counts_I151', 'counts_I348', 'counts_I701']
ICD_names += ['counts_I700', 'counts_I311', 'counts_R71']
ICD_names += ['counts_I099', 'counts_I770', 'counts_I051']
ICD_names += ['counts_I250', 'counts_I400', 'counts_I238']
ICD_names += ['counts_I484', 'counts_I513', 'counts_I510']
ICD_names += ['counts_I708', 'counts_I742', 'counts_I653']
ICD_names += ['counts_I721', 'counts_I058', 'counts_I159']
ICD_names += ['counts_I729', 'counts_I820', 'counts_I241']
ICD_names += ['counts_I310', 'counts_I722', 'counts_I7090']
ICD_names += ['counts_I741', 'counts_I772', 'counts_I716']
ICD_names += ['counts_I516', 'counts_I288', 'counts_I669']
ICD_names += ['counts_I515', 'counts_I7010', 'counts_I060']
ICD_names += ['counts_I089', 'counts_I740', 'counts_I711', 'counts_I139']
ICD_names = np.array(ICD_names)

real_names = ["hypertension", "chest pain", "haematuria"]
real_names += ["nausea/vomiting", "unspecified morbidity", "syncope and collapse"]
real_names += ["retention of urine", "atherosclerotic HD", "atrial fibrillation"]
real_names += ["haemorrhoids", "headache", "varicose veins"]
real_names += ["pleural effusion", "chest pain", "dizziness and giddiness"]
real_names += ["pulmonary embolism", "hypotension", "atrial fibrillation"]
real_names += ["angina", "chronic ischaemic HD", "any_HF"]
real_names += ["dysphagia", "LV_HF", "myocardial infarction"]
real_names += ["polyuria", "palpitations", "dyspnoea"]
real_names += ["urinary incontinence", "malaise and fatigue", "phlebitis"]
real_names += ["cough", "orthostatic hypotension", "bradycardia"]
real_names += ["haemorrhoids", "faecal incontinence", "atrial fibrillation"]
real_names += ["haemorrhoids", "myocardial infarction", "abnormal lung imaging"]
real_names += ["chronic ischaemic HD", "myocardial infarction", "cerebral infarction"]
real_names += ["angina", "peripheral vascular disease", "precordial pain"]
real_names += ["haemoptysis", "Unknown_HF", "nonspecific rash"]
real_names += ["haemorrhoids", "tachycardia", "tachycardia"]
real_names += ["cong_HF'", "valve stenosis", "raynaud's syndrome", "CMyo"]
real_names += ["epistaxis", "flatulence", "ascites"]
real_names += ["heartburn", "stroke", "stroke"]
real_names += ["cardiomegaly", "myocardial infarction", "myocardial infarction"]
real_names += ["mouth breathing", "cerebrovascular disease", "heart block"]
real_names += ["senility", "heart block", "Cardiac murmur"]
real_names += ["cerebrovascular disease", "brain haemorrhage", "varicose veins"]
real_names += ["high blood pressure reading", "thrombosis", "valve insufficiency"]
real_names += ["heart disease", "cerebral infarction", "ischaemic heart disease"]
real_names += ["haemorrhoids", "venous insufficiency", "atherosclerosis"]
real_names += ["jaundice", "mitral valve prolapse", "cardiac arrhythmia"]
real_names += ["gangrene", "Stricture of artery", "oesophageal varices"]
real_names += ["tachycardia", "cerebral infarction", "phlebitis"]
real_names += ["hypertensive renal disease", "low blood pressure reading", "valve stenosis"]
real_names += ["heart block", "cardiac arrest", "oedema"]
real_names += ["phlebitis", "breathing abnormalities", "cerebral aneurysm"]
real_names += ["brain haemorrhage", "angina", "occlusion and stenosis"]
real_names += ["valve disease", "heart block", "lymphoedema"]
real_names += ["atherosclerosis", "myocardial infarction", "naevus"]
real_names += ["oesophageal varices", "ischaemic heart disease", "atrial fibrillation"]
real_names += ["varicose veins", "phlebitis", "cardiac depolarisation"]
real_names += ["cardiac arrest", "pain in throat", "valve disease"]
real_names += ["pericardial effusion", "aortic aneurysm", "ventricular fibrillation"]
real_names += ["cardiac depolarisation", "pleurisy", "preexcitation syndrome"]
real_names += ["oesophageal varices", "valve insufficiency", "pleural effusion"]
real_names += ["proteinuria", "splenomegaly", "cerebral infarction"]
real_names += ["wheezing", "pulmonary hypertension", "heart block"]
real_names += ["valve insufficiency", "endocarditis", "brain haemorrhage"]
real_names += ["atherosclerosis", "disorder of arteries", "hypertensive heart disease"]
real_names += ["pulmonary heart disease", "heart block", "myocardial infarction"]
real_names += ["aortic aneurysm", "pulmonary oedema", "cerebrovascular disease"]
real_names += ["valve disease", "abnormal sputum", "hyperventilation"]
real_names += ["myocardial infarction", "arteritis", "sick sinus syndrome"]
real_names += ["atherosclerosis", "aortic aneurysm", "heart block"]
real_names += ["brain haemorrhage", "brain haemorrhage", "myocardial infarction"]
real_names += ["phlebitis", "cerebral infarction", "pulmonary hypertension"]
real_names += ["valve disease", "varices", "valve disease"]
real_names += ["chest pain on breathing", "thrombosis", "heart disease"]
real_names += ["cardiorespiratory symptoms", "abnormal white blood cells", "pericardium disease"]
real_names += ["thrombosis", "atrial fibrillation", "heart block"]
real_names += ["valve disease", "oedema", "tachycardia"]
real_names += ["cerebral infarction", "brain haemorrhage", "heart beat abnormalities"]
real_names += ["cerebral infarction", "heart block", "arrhythmia"]
real_names += ["lymphangitis", "valve disease", "atrial flutter"]
real_names += ["thrombosis ", "arrhythmia", "varices"]
real_names += ["respiratory arrest", "heart block", "hepatomegaly"]
real_names += ["valve disease", "brain haemorrhage", "brain haemorrhage"]
real_names += ["vein disease", "heart block", "pericarditis"]
real_names += ["thrombosis", "hypotension (drugs)", "lymphatic disease"]
real_names += ["angina with spasm", "brain haemorrhage", "atherosclerosis"]
real_names += ["aneurysm", "varices", "asphyxia"]
real_names += ["endocarditis", "angiopathy", "stridor"]
real_names += ["oesophageal varices", "circulatory disease", "valve insufficiency"]
real_names += ["valve stenosis", "hypotension", "renal disease"]
real_names += ["pulmonary embolism", "valve stenosis", "haemorrhagic telangiectasia"]
real_names += ["chordae tendineae rupture", "thrombosis", "cardiac depolarisation"]
real_names += ["hepatomegaly and splenomegaly", "aneurysm", "cachexia"]
real_names += ["brain haemorrhage", "compression of vein", "thrombosis"]
real_names += ["cerebral infarction", "heart block", "aortic dissection"]
real_names += ["aortic aneurysm", "haemorrhage", "thrombosis"]
real_names += ["capillary disease", "atherosclerosis", "hypotension"]
real_names += ["heart block", "myocardial infarction", "anuria and oliguria"]
real_names += ["valve stenosis", "disorder of arteries", "aortic aneurysm"]
real_names += ["postprocedural disorders", "valve disease", "valve insufficiency"]
real_names += ["atherosclerosis", "phlebitis", "myocarditis"]
real_names += ["pulmonary heart disease", "thrombosis", "haemorrhage"]
real_names += ["aneurysm", "urethral discharge", "heart block"]
real_names += ["thrombosis", "thrombosis", "postprocedural disorders"]
real_names += ["artery stenosis", "oedema", "aneurysm"]
real_names += ["cerebrovascular disease", "peripheral vascular disease", "myocardial infarction"]
real_names += ["secondary hypertension", "valve disease", "atherosclerosis"]
real_names += ["atherosclerosis", "pericarditis", "abnormal red blood cells"]
real_names += ["rheumatic illness", "arteriovenous fistula", "rheumatic illness"]
real_names += ["atherosclerosis", "myocarditis (infection)", "myocardial infarction"]
real_names += ["atrial flutter", "thrombosis", "septal defect"]
real_names += ["Atherosclerosis", "thrombosis", "occlusion and stenosis"]
real_names += ["aneurysm", "valve disease", "secondary hypertension"]
real_names += ["aneurysm", "budd-chiari syndrome", "dressler's syndrome"]
real_names += ["pericarditis", "aneurysm", "atherosclerosis"]
real_names += ["thrombosis", "rupture of artery", "aortic aneurysm"]
real_names += ["heart disease", "pulmonary vessel disease", "occlusion and stenosis"]
real_names += ["myocardial degeneration", "atherosclerosis", "rheumatic illness"]
real_names += ["valve disease", "thrombosis", "aortic aneurysm", "hypertensive heart disease"]
real_names = np.array(real_names)
real_names_new_inds = np.array([np.where(ICD_names == name)[0][0] for name in sorted_Y_names])
if np.all(ICD_names[real_names_new_inds] == sorted_Y_names):
    ICD_names = ICD_names[real_names_new_inds]
    real_names = real_names[real_names_new_inds]
else:
    print("exiting: sorting error on line 365")
    exit()

# https://www.ncbi.nlm.nih.gov/books/NBK535419/ https://www.sciencedirect.com/science/article/pii/S0735109720377755
# Heart disease (HD): severe heart conditions: everything from an atherosclerotic heart attack to left ventricular heart failure. 
# artery and vein disease (AVD): diseases the blood vessels: everything from atherosclerosis to haemorrhoids
# abnormal heart activity (AHA): what it sounds like: high/low blood pressure, atrial fibrillation, etc  
# possible outcomes of cardiovascular disease (POC): anything that is sometimes caused by CHD, PAD, and AHA: everything from coughing to strokes to urine retention
# misc: anything rarely caused by heart disease: Dizziness and giddiness and Haemoptysis (coughing blood) would be examples

categories = ["AHA", "POC", "MISC"]
categories += ["MISC", "MISC", "POC"]
categories += ["POC", "HD", "AHA"]
categories += ["AVD", "POC", "AVD"]
categories += ["POC", "POC", "MISC"]
categories += ["AVD", "AHA", "AHA"]
categories += ["POC", "HD", "HD"]
categories += ["MISC", "HD", "HD"]
categories += ["POC", "AHA", "POC"]
categories += ["MISC", "POC", "AVD"]
categories += ["POC", "AHA", "AHA"]
categories += ["AVD", "MISC", "AHA"]
categories += ["AVD", "HD", "POC"]
categories += ["HD", "HD", "HD"]
categories += ["POC", "AVD", "POC"]
categories += ["MISC", "HD", "MISC"]
categories += ["AVD", "AHA", "AHA"]
categories += ["HD", "AVD", "AVD", "HD"]
category_names = real_names[:55]

# TODO: Make a table with the elements of each category:
names = np.array(["POC", "HD", "AVD", "AHA", "MISC"]).reshape(-1, 1)
counts = np.sum(np.array(categories) == names, axis = 1).astype(float)
percentages = (counts/np.sum(counts))*100
labels = [n[0] + " (" + str(np.round(p, 1)) + "%)" for n, p in zip(names, percentages)]
colors_vec = [[0, 1, 1], [0, 0.80, 1], [0, 1, 1], [0, 0.8, 1], [0.8, 0, 0.3]]
plt.pie(counts, labels = labels, labeldistance = 0.25, colors = colors_vec, rotatelabels = True, counterclock = False, startangle = 90)
plt.savefig("aaa.png")
plt.clf()

ICD_to_real = dict(zip(ICD_names, real_names))
N = 6
if not os.path.isdir("pie_charts"):
    os.mkdir("pie_charts")
for i in range(k):
    
    Y_roti = COPY(Y_rots[i, :])

    sorted_indices = np.flip(np.argsort(Y_roti))
    sorted_Y_roti = Y_roti[sorted_indices]
    sorted_Y_roti_standardized = sorted_Y_roti/np.sum(np.abs(sorted_Y_roti))
    sorted_ICD_names = y_cols[1:][sorted_indices]
    sorted_real_names = np.array([ICD_to_real[name] for name in sorted_ICD_names])

    N1 = np.min([N, np.sum(sorted_Y_roti_standardized > 0)])
    N2 = np.min([N, np.sum(sorted_Y_roti_standardized < 0)])
    pie_labs_pos = sorted_real_names[:N1]
    pie_labs_neg = sorted_real_names[-N2:]
    sorted_Y_roti_pos = sorted_Y_roti_standardized[:N1]
    sorted_Y_roti_neg = sorted_Y_roti_standardized[-N2:]
    pie_vals_pos = sorted_Y_roti_pos/np.sum(sorted_Y_roti_pos)
    pie_vals_neg = sorted_Y_roti_neg/np.sum(sorted_Y_roti_neg)
    col_vec_neg = np.array(2*[[0.8, 0, 0], [0.8, 0, 0.5], [0.8, 0, 0.7]])
    col_vec_pos = np.array(2*[[0, 0.7, 1], [0, 0.85, 1], [0, 1, 1]])
    pie_labs_pos = sorted_real_names[:N1]
    pie_labs_neg = sorted_real_names[-N2:]
    if i == 1:
        sorted_real_names[-N2:-1] = "other"

    pos_df = pd.DataFrame(np.array([pie_labs_pos, pie_vals_pos], dtype = object).T).groupby(0, as_index=False).sum()
    neg_df = pd.DataFrame(np.array([pie_labs_neg, pie_vals_neg], dtype = object).T).groupby(0, as_index=False).sum()

    # figsize: width, length
    if N2 > 0: 
        fig, ax = plt.subplots(1,2, figsize=(30, 15))
        ax[0].pie(pos_df[1], labels = pos_df[0], radius = 1.45, textprops={'fontsize': 24}, labeldistance = 0.20, colors = col_vec_pos[:len(pos_df)], rotatelabels = True, counterclock = False, startangle = 90)
        ax[1].pie(neg_df[1], labels = neg_df[0], radius = 1.45, textprops={'fontsize': 24}, labeldistance = 0.20, colors = col_vec_neg[:len(neg_df)], rotatelabels = True, counterclock = False, startangle = 90)
        pos_percent = np.sum(np.sum(sorted_Y_roti[:N]))/(np.sum(sorted_Y_roti[:N]) + np.sum(np.abs(sorted_Y_roti[-N:])))
        neg_percent = np.sum(np.sum(np.abs(sorted_Y_roti[-N:])))/(np.sum(sorted_Y_roti[:N]) + np.sum(np.abs(sorted_Y_roti[-N:])))
        title0 = 'top ' + str(N1) + ' positive contributions (' + str(np.round(pos_percent*100, 1))[:4] + '%)\nto latent phenotype ' + str(i + 1) 
        title1 = 'top ' + str(N2) + ' negative contributions (' + str(np.round(neg_percent*100, 1))[:4] + '%)\nto latent phenotype ' + str(i + 1)
        ax[0].set_title(title0, fontsize = 40, y = ax[0].title._y + 0.1)
        ax[1].set_title(title1, fontsize = 40, y = ax[1].title._y + 0.1)
        fig.subplots_adjust(bottom=0.05, top=0.85)
        fig.savefig("pie_charts/feature" + str(i + 1) + ".png")
    else:
        fig, ax = plt.subplots(1,1, figsize=(30, 15))
        ax.pie(pos_df[1], labels = pos_df[0], radius = 1.30, textprops={'fontsize': 24}, labeldistance = 0.20, colors = col_vec_pos[:len(pos_df)], rotatelabels = True, counterclock = False, startangle = 90)
        title0 = 'top ' + str(N1) + ' contributions \nto latent phenotype ' + str(i + 1) 
        ax.set_title(title0, fontsize = 40, y = ax.title._y + 0.05)
        fig.subplots_adjust(bottom=0.05, top=0.85)
        fig.savefig("pie_charts/feature" + str(i + 1) + ".png")


'''
path = "logistic_SVD_output/phenotypes" + str(k) + ".txt"
phenotypes = pd.read_csv(path, delimiter = "\t").to_numpy()[:, 1:]
r_vals, p_vals = [], []
for x in COPY(X.T):
    x[np.isnan(x)] = np.nanmean(x)
    path = "logistic_SVD_output/phenotypes" + str(k) + ".txt"
    phenotypes = pd.read_csv(path, delimiter = "\t").to_numpy()[:, 1:]
    corrs = [pearsonr(x, p) for p in phenotypes.T]
    r_vals.append([i[0] for i in corrs])
    p_vals.append([i[1] for i in corrs])

pdb.set_trace()
r_vals, p_vals = np.array(r_vals), np.array(p_vals)
p_bon = 0.05/len(X[0])
p_vals_for_pos_R = p_vals[r_vals[:, 1] > 0, 1]
num_p_vals_for_pos_R = np.sum(p_vals_for_pos_R < p_bon)
p_vals_for_neg_R = p_vals[r_vals[:, 1] < 0, 1]
num_p_vals_for_neg_R = np.sum(p_vals_for_neg_R < p_bon)

top_pos_terms = X_cols[np.argsort(r_vals[:, 1])[-num_p_vals_for_pos_R:]]
top_neg_terms = X_cols[np.argsort(r_vals[:, 1])[:num_p_vals_for_neg_R]]
# 6144-4.0: never eat added sugar 23105-0.0: metabolic rate 894-average: exercise 6177-1.0: BP, chol meds
neg_weird = np.array(["pack-years", "6177-2.0", "6177-1.0'", "894-average", "23105-0.0", "6144-4.0", "annual-consumption", "21001-0.0"])
# not a victim of violence, been in confiding relationship, not assaulted, never witnessed violent death, 
neg_expected = np.array(["20529-0.0", "20522-average", "20531-0.0", "20530-0.0"])
# very bad relationship, very bad childhood, bad childhood 2, bad childhood 3, bad relationship 2
pos_expected = np.array(["20524-average", "20490-average", "20488-average", "20487-average", "20523-average"])
expected = np.concatenate([neg_expected, pos_expected])
X2 = COPY(X.T)
for x in X2: x[np.isnan(x)] = np.nanmean(x)
X2 = X2.T
X_weird = np.concatenate([X2[:, np.isin(X_cols, neg_weird)], np.ones((len(X2), 1))], axis = 1)
X_expected = np.concatenate([X2[:, np.isin(X_cols, expected)], np.ones((len(X2), 1))], axis = 1)
model_weird = sm.OLS(phenotypes[:, 1], X_weird)
model_expected = sm.OLS(phenotypes[:, 1], X_expected)
model_full = sm.OLS(phenotypes[:, 1], X2)
weird_results = model_weird.fit()
expected_results = model_expected.fit()
expected_full = model_full.fit()
#model = sm.OLS(phenotype, X)
#model_sub = sm.OLS(phenotype, X_sub)
#model_results = model.fit()
#model_sub_results = model_sub.fit()
#p_values.append(model_results.compare_lr_test(model_sub_results)[1])
'''