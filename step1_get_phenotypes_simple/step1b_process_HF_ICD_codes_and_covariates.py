import numpy as np
import pandas as pd
from copy import deepcopy as COPY
from scipy.stats import pearsonr
from tqdm import tqdm
import pdb

cong_HF = ["counts_" + ICD for ICD in ["I110", "I130", "I132", "I500"]]
CMyo = ["counts_" + ICD for ICD in ["I255", "I420", "I425", "I426", "I427", "I428", "I429", "I431"]]
LV_HF = "counts_I501"
Unknown_HF = "counts_I509"
file1 = pd.read_csv("40001for_HF.txt", delimiter = "\t", header = 0).astype(int)
file2 = pd.read_csv("40002for_HF.txt", delimiter = "\t", header = 0).astype(int)
file3 = pd.read_csv("41270for_HF.txt", delimiter = "\t", header = 0).astype(int)
case1 = np.all(file1["eid"] == file2["eid"])
case2 = np.all(file2["eid"] == file3["eid"])
if not np.logical_and(case1, case2):
    print("error: the eids from file1, file2, and file3 should be sorted identically")
    exit()
y_data = COPY(file1.loc[:, :])
y_data.loc[:, y_data.columns != "eid"] += file2.loc[:, file2.columns != "eid"]
y_data.loc[:, y_data.columns != "eid"] += file3.loc[:, file3.columns != "eid"]
y_data.loc[:, y_data.columns != "eid"] = y_data.loc[:, y_data.columns != "eid"] > 0
y_data["cong_HF"] = np.any(y_data.loc[:, np.isin(y_data.columns, cong_HF)].to_numpy(), axis = 1)
y_data["CMyo"] = np.any(y_data.loc[:, np.isin(y_data.columns, CMyo)].to_numpy(), axis = 1)
y_data["LV_HF"] = np.any(y_data.loc[:, np.isin(y_data.columns, LV_HF)].to_numpy(), axis = 1)
y_data["Unknown_HF"] = np.any(y_data.loc[:, np.isin(y_data.columns, Unknown_HF)].to_numpy(), axis = 1)
all_HF = ["cong_HF", "CMyo", "LV_HF", "Unknown_HF"]
y_data["any_HF"] = np.any(y_data.loc[:, np.isin(y_data.columns, all_HF)].to_numpy(), axis = 1)

HF = ['cong_HF', 'CMyo', 'LV_HF', 'Unknown_HF', 'any_HF']
y_data = y_data.loc[:, ["eid"] + HF]

special_feature_path = "/home/greggj/pleiotropy_and_GxE/step1_get_my_phenotypes_complex/reformatted_fields.txt"
special_fields = ["eid", "max-46-47", "is_urban", "annual-consumption", "pack-years", "102-average", "4079-average", "4080-average"]
# duration of walks, moderate activity, vigorous activity
special_fields += ["874-average", "894-average", "914-average"]
special_fields += ["20487-average", "20488-average", "20490-average", "20491-average", "20521-average", "20522-average", "20523-average"]
special_fields += ["20524-average", "20525-average", "20526-0.0", "20526-1.0", "20526-2.0"]
special_fields += ["20527-0.0", "20527-1.0", "20527-2.0", "20529-0.0", "20529-1.0", "20529-2.0"]
special_fields += ["20530-0.0", "20530-1.0", "20530-2.0", "20531-0.0", "20531-1.0", "20531-2.0"]
special_fields += ["6160-1.0", "6160-2.0", "6160-3.0", "6160-4.0", "6160-5.0"]
special_fields += ["6177-1.0", "6177-2.0", "6177-3.0", "6177-4.0", "6177-5.0"]
special_fields += ["6179-1.0", "6179-2.0", "6179-3.0", "6179-4.0", "6179-5.0", "6179-6.0"]
special_fields += ["6155-1.0", "6155-2.0", "6155-3.0", "6155-4.0", "6155-5.0", "6155-6.0", "6155-7.0"]
special_fields += ["6145-1.0", "6145-2.0", "6145-3.0", "6145-4.0", "6145-5.0", "6145-6.0"]
special_fields += ["6144-1.0", "6144-2.0", "6144-3.0", "6144-4.0", "6144-5.0"]
special_fields += ["6138-1.0", "6138-2.0", "6138-3.0", "6138-4.0", "6138-5.0", "6138-6.0"]
special_fields += ["1767-0.0", "1389-average", "1438-average", "1458-average", "1478-average", "1488-average"]
special_fields += ["1498-average", "1528-average", "1349-average", "1408-average", "1359-average"]
special_fields += ["1369-average", "1379-average", "1319-average", "1329-average", "1339-average"]
special_fields += ["1289-average", "1299-average", "1309-average", "1110-average", "699-average"]
special_fields += ["670-1.0", "670-2.0", "670-3.0", "670-4.0", "670-5.0"]
special_fields += ["680-1.0", "680-2.0", "680-3.0", "680-4.0", "680-5.0", "680-6.0", "189-average"]
special_fields += ["1687-1.0", "1687-2.0", "1687-3.0", "1697-1.0",  "1697-2.0",  "1697-3.0"]
special_fields += ["2306-0.0", "2306-2.0", "2306-3.0", "1031-average", "24506-average"]

special_features = pd.read_csv(special_feature_path, delimiter = "\t", usecols = special_fields, header = 0)

X_eids = pd.read_csv("ICD_eids.txt", delimiter = "\t", header = 0)
X_ICD_main = pd.read_csv("ICD_data_main.txt", delimiter = "\t", header = 0)
X_colnames = X_ICD_main.columns
X_ICD_other_PCs = pd.read_csv("ICD_data_other_PCs.txt", delimiter = "\t", header = None, usecols = np.arange(50))
X_ICD_other_PCs.columns = ["ICD_PC" + str(i) for i in X_ICD_other_PCs.columns]
X_ICD_dates_main = pd.read_csv("ICD_data_dates_main.txt", delimiter = "\t", header = 0)
X_ages_colnames = X_ICD_dates_main.columns
X_ICD_dates_other_PCs = pd.read_csv("ICD_data_dates_other_PCs.txt", delimiter = "\t", header = None, usecols = np.arange(50))
X_ICD_dates_other_PCs.columns = ["ICD_date_PC" + str(i) for i in X_ICD_dates_other_PCs.columns]
X_data = pd.concat([X_eids, 
                    X_ICD_main.reindex(X_eids.index), 
                    X_ICD_other_PCs.reindex(X_eids.index), 
                    X_ICD_dates_main.reindex(X_eids.index), 
                    X_ICD_dates_other_PCs.reindex(X_eids.index)], axis=1)

other_fields = ["1920", "1930", "1940", "1950", "1960", "1970", "1980"]
other_fields += ["1990", "2000", "2010", "2020", "2030", "2040", "2050", "2060"]
other_fields += ["2335", "2443", "21001", "22001", "24508", "23105"]
other_fields += ["30600", "30610", "30620", "30630", "30640", "30650", "30660", "30670", "30680"]
other_fields += ["30690", "30700", "30710", "30720", "30730", "30740", "30750", "30760", "30770"]
other_fields += ["30780", "30790", "30810", "30830", "30840", "30850", "30860", "30870", "30880", "30890"]
other_fields += ["24003", "24004", "24005", "24006", "24007", "24008", "24009"]
other_fields += ["24010", "24011", "24012", "24015", "24018", "24019", "24020", "24021", "24022"]
other_fields += ["30000", "30010", "30020", "30030", "30040", "30050", "30060", "30070", "30080"] 
other_fields += ["30090", "30100", "30110", "30120", "30130", "30140", "30150", "30160", "30170"]
other_fields += ["30180", "30190", "30200", "30210", "30220", "30230", "30240"]
other_fields += ["30250", "30260", "30270", "30280", "30290", "30300"]
other_fields = [name + "-0.0" for name in other_fields]
other_fields += special_fields

HF_metadata = pd.read_csv("HF_metadata.txt", delimiter = "\t", header = 0)
HF_metadata = HF_metadata.merge(special_features, how = "inner", on = "eid")
metadata_cols = ["eid", "heart attack", "angina", "stroke", "high BP", "age", "22001-0.0"] + other_fields
X_metadata = HF_metadata[metadata_cols]
X_metadata_colnames = X_metadata.columns

# removes features with too little variance
all_data = X_data.merge(HF_metadata, on = "eid", how = "inner").astype(float)
good_col_indices = (np.sum(np.logical_and(all_data != 0, np.isnan(all_data) == False), axis = 0) > 30).to_numpy()
good_col_indices = np.logical_and(good_col_indices, np.nanvar(all_data, axis = 0) > 0)
good_cols = all_data.columns[good_col_indices]
all_data = all_data[good_cols]

# Arguable fields removed:
# cholesteral and other biomarkers: they tend to be highly heritable. Non-heritable portions should be contained in other covariates. 
# field number order: risk taking, BMI, gender, basal metabolic rate
# 24003-24022: various polution indices
# distance to coast, 
environmental_cols = ['eid', 'private insurance degree', 'private insurance','C1', 'C2', 
                      '2040-0.0', '21001-0.0', '22001-0.0', '23105-0.0',  
                      '24003-0.0', '24004-0.0', '24005-0.0', '24006-0.0', 
                      '24007-0.0', '24008-0.0', '24009-0.0', '24010-0.0', 
                      '24011-0.0', '24012-0.0', '24015-0.0', '24018-0.0', 
                      '24019-0.0', '24020-0.0', '24021-0.0', '24022-0.0', 
                      '24508-0.0', 'age', 'income', 'pack-years', 'is_urban', 'annual-consumption',
                      "874-average", "894-average", "914-average", "20487-average", "20488-average", 
                      "20490-average", "20491-average", "20521-average", "20522-average", 
                      "20523-average", "20524-average", "20525-average", "20526-0.0", "20526-1.0", "20526-2.0",
                      "20527-0.0", "20527-1.0", "20527-2.0", "20529-0.0", "20529-1.0", "20529-2.0",
                      "20530-0.0", "20530-1.0", "20530-2.0", "20531-0.0", "20531-1.0", "20531-2.0",
                      "6160-1.0", "6160-2.0", "6160-3.0", "6160-4.0", "6160-5.0", 
                      "6177-1.0", "6177-2.0", "6177-3.0", "6177-4.0", "6177-5.0",
                      "6179-1.0", "6179-2.0", "6179-3.0", "6179-4.0", "6179-5.0", "6179-6.0",
                      "6155-1.0", "6155-2.0", "6155-3.0", "6155-4.0", "6155-5.0", "6155-6.0", "6155-7.0",
                      "6145-1.0", "6145-2.0", "6145-3.0", "6145-4.0", "6145-5.0", "6145-6.0",
                      "6144-1.0", "6144-2.0", "6144-3.0", "6144-4.0", "6144-5.0",
                      "6138-1.0", "6138-2.0", "6138-3.0", "6138-4.0", "6138-5.0", "6138-6.0",]
special_fields += ["1767-0.0", "1389-average", "1438-average", "1458-average", "1478-average", "1488-average"]
special_fields += ["1498-average", "1528-average", "1349-average", "1408-average", "1359-average"]
special_fields += ["1369-average", "1379-average", "1319-average", "1329-average", "1339-average"]
special_fields += ["1289-average", "1299-average", "1309-average", "1110-average", "699-average"]
special_fields += ["670-1.0", "670-2.0", "670-3.0", "670-4.0", "670-5.0"]
special_fields += ["680-1.0", "680-2.0", "680-3.0", "680-4.0", "680-5.0", "680-6.0", "189-average"]
special_fields += ["1687-average", "1697-average", "2306-average", "1031-average", "24506-average"]


all_data = all_data[environmental_cols]

# no features extremely correlated with heart failure were found.
HF_correlations = []
for ycol in y_data[HF].to_numpy(dtype = float).T:
    correlation_set = []
    X = all_data[all_data.columns[all_data.columns != "eid"]].to_numpy(dtype = float).T
    i = -1
    for xcol in tqdm(X):
        i += 1
        val_indices = np.isnan(xcol) == False
        r = pearsonr(xcol[val_indices], ycol[val_indices])[0]
        correlation_set.append(pearsonr(xcol[val_indices], ycol[val_indices])[0])
    HF_correlations.append(correlation_set)
HF_correlations = np.array(HF_correlations)

all_data.to_csv("X.txt", sep = "\t", header = True, index = False)
y_data.to_csv("y.txt", sep = "\t", header = True, index = False)


