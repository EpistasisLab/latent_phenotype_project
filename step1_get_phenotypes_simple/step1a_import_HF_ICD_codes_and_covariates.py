import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
from copy import deepcopy as COPY
from functools import reduce
from tqdm import tqdm
import time
import pdb
import os

def is_field(col_name, fields):

    """
    Purpose
    -------
    Parameters
    ----------
    Returns
    -------
    """

    status = False
    column_has_0th_instance = False

    # These lines assume standard colname notation of "field-X.Y", where X is the instance number, and Y is the rep number.
    partial_col_names = [field + "-" for field in fields]
    colname_has_partial = np.any([col_name[0:(len(part))] == part for part in partial_col_names])
    if colname_has_partial: column_has_0th_instance = col_name.split("-")[1][0] == "0"

    # Also, if there is only 1 instance and 1 rep, then the colname is just "field"
    if column_has_0th_instance or col_name in fields: status = True
    
    return(status)

# 22006, 21000 = is white british, self-declared ethnicity
# 40001, 40002 = primary, secondary cause of death
# 41270, 41280 = unique ICD codes, date of first ICD code instance
# 34, 52 = birth year, birth month
# 4674, 738 = private insurance, household income
# 41249, 41245, 41246 = reason for hosptial visit, attending doctor's specialty, specialty being practiced
# 21003, 6150 age when attending center, misc heart problems diagnosed by doctor
# 129 (north) or 130 (east) coordinate, 189 townsend deprevation index 
# 699 = length of time at current address
# 864, 884, 904 = exercise frequencies
# 1200, 1160 = subjectives, objective sleeplessness
# 1920 - 2060 in intervals of 10: negative mental traits 
# 2335, 2443, 21001, 22001 = chest pain, diabetes not ICD codes, BMI, gender
# 24508, 23105 = distance to coast, basal metabolic rate
#30600 - 30890 in intervals of 10, biomarkers
#24003 - 24022 with some missing: various polution measures. 
#30000 - 30300 in intervals of 10, cell markers
my_fields = ["eid", "20002", "22006", "21000", "40002", "40001", "41270", "41280", "34", "52"]
my_fields += ["4674", "738", "41249", "41245", "41246", "21003", "6150"]
my_fields += ["129", "130", "699", "864", "884", "904", "1200", "1160"]
my_fields += ["1920", "1930", "1940", "1950", "1960", "1970", "1980", "1990", "2000", "2010", "2020", "2030", "2040", "2050", "2060"]
my_fields += ["2335", "2443", "21001", "22001", "24508", "23105"]
my_fields += ["30600", "30610", "30620", "30630", "30640", "30650", "30660", "30670", "30680"]
my_fields += ["30690", "30700", "30710", "30720", "30730", "30740", "30750", "30760", "30770"]
my_fields += ["30780", "30790", "30810", "30830", "30840", "30850", "30860", "30870", "30880", "30890"]
my_fields += ["24003", "24004", "24005", "24006", "24007", "24008", "24009"]
my_fields += ["24010", "24011", "24012", "24015", "24018", "24019", "24020", "24021", "24022"]
my_fields += ["30000", "30010", "30020", "30030", "30040", "30050", "30060", "30070", "30080"] 
my_fields += ["30090", "30100", "30110", "30120", "30130", "30140", "30150", "30160", "30170"]
my_fields += ["30180", "30190", "30200", "30210", "30220", "30230", "30240"]
my_fields += ["30250", "30260", "30270", "30280", "30290", "30300"]

if os.path.isfile("phenotype_info.txt") == False:
    #path1 = "/project/UKB_moore/UKB_50978/phenotype/penn_freeze_11132019/ukb38304.csv"
    path1 = "/project/UKB_moore/UKB_50978/phenotype/penn_freeze_05142021/ukb46981.csv"
    path2 = "/project/UKB_moore/UKB_50978/phenotype/penn_add_08042020/ukb43023.csv"
    path3 = "/project/UKB_moore/UKB_50978/phenotype/penn_add_12082020/ukb44783.csv"
    path4 = "/project/UKB_moore/UKB_50978/phenotype/penn_add_12232020/ukb44910.csv"
   
    fields1 = pd.read_csv(path1, delimiter = ',', usecols = lambda col_name: is_field(col_name, my_fields), header = 0, dtype = str)
    fields2 = pd.read_csv(path2, delimiter = ',', usecols = lambda col_name: is_field(col_name, my_fields), header = 0, dtype = str)
    fields3 = pd.read_csv(path3, delimiter = ',', usecols = lambda col_name: is_field(col_name, my_fields), header = 0, dtype = str)
    fields4 = pd.read_csv(path4, delimiter = ',', usecols = lambda col_name: is_field(col_name, my_fields), header = 0, dtype = str)

    fields = fields1.merge(fields2, how = "inner", on = "eid")
    fields = fields.merge(fields3, how = "inner", on = "eid")
    fields = fields.merge(fields4, how = "inner", on = "eid")
    fields.to_csv("phenotype_info.txt", sep = "\t", header = True, index = False)
else:
    fields = pd.read_csv("phenotype_info.txt", delimiter = "\t", header = 0,  dtype = str)

is_white = np.isin(fields["21000-0.0"], ["1", "1001", "1002", "1003"])
self_declared_illness_cols = fields.columns[np.array(["20002" in col for col in fields.columns])]
has_HCM1 = np.any(fields[self_declared_illness_cols] == "1588" , axis = 1)
has_HCM2 = np.any(np.isin(fields, ["I421", "I422"]), axis = 1)
has_no_HCM = np.logical_or(has_HCM1, has_HCM2) == False
fields = fields.loc[np.logical_and(is_white, has_no_HCM), :]
fields[["eid", "eid"]].to_csv("../step2_get_UKB_samples/eids.tab", sep = "\t", header = False, index = False)

#---------------------------------------------------------------------------------------------------------------------------------
# "fields" has all non-white and hypertrophic CM individuals removed
# At the bottom, I carve out a dataset for lstm dimensionality reduction
#---------------------------------------------------------------------------------------------------------------------------------

HF_ICD_codes = ["I110", "I130", "I132", "I255", "I420", "I425", "I428", "I429", "I500", "I501", "I509"]
colnames = fields.columns
colfields = np.array([name.split("-")[0] for name in colnames])
ICD_codes_1 = fields.loc[:, colnames[np.isin(colfields, ["eid", "40002"])]]
ICD_codes_2 = fields.loc[:, colnames[np.isin(colfields, ["eid", "40001"])]]
ICD_codes_3 = fields.loc[:, colnames[np.isin(colfields, ["eid", "41270"])]]
ICD_codes_all = ICD_codes_1.merge(ICD_codes_2, how = "inner", on = "eid")
ICD_codes_all = ICD_codes_all.merge(ICD_codes_3, how = "inner", on = "eid")
ICD_codes_all[ICD_codes_all.isna()] = "NA"
del ICD_codes_all["eid"]
ICD_codes_all = ICD_codes_all.to_numpy(dtype = str).reshape(-1)
ICD_codes_all = np.unique(ICD_codes_all[ICD_codes_all != "NA"])

ICD_codes = np.setdiff1d(ICD_codes_all, HF_ICD_codes)
ICD_codes = np.sort(ICD_codes)
main_indices1 = [code[0] == "I" for code in ICD_codes]
main_indices2 = [code[0] == "R" and int(code[1:4]) <= 99 for code in ICD_codes]
main_ICD_codes = ICD_codes[np.logical_or(main_indices1, main_indices2)].tolist()
main_ICD_codes += ['J81', 'J90', 'J91', 'R160', 'R161', 'R162', 'R600', 'R601', 'R609']
main_ICD_codes = np.sort(main_ICD_codes)

birthdates = pd.to_datetime(fields.loc[:, "34-0.0"] + "-" + fields.loc[:, "52-0.0"] + "-15")
for col in colnames[colfields == "41280"]:
    sickdates = pd.to_datetime(fields.loc[:, col])
    fields[col] = (sickdates - birthdates).dt.total_seconds()/(24*3600*365)

ICD_counts_subsets = []
HF_ICD_counts_subsets = []
ICD_dates_original = fields.loc[:, colnames[np.isin(colfields, ["eid", "41280"])]].to_numpy()
for field in ["40002", "40001", "41270"]:
    subset = fields.loc[:, colnames[np.isin(colfields, ["eid", field])]]
    subset[subset.isna()] = "NA"
    subset_colnames = subset.columns
    subset_np = subset.to_numpy()
    ICD_counts = np.zeros((len(subset), (len(ICD_codes) + 1)))
    ICD_dates = np.zeros((len(subset), (len(ICD_codes) + 1)))
    HF_ICD_counts = subset.loc[:, ["eid"]]
    HF_ICD_dates = subset.loc[:, ["eid"]]
    for i in tqdm(range(len(subset_np))):
        row = subset_np[i]
        row = row[row != "NA"]
        sorted_row_indices = np.argsort(row[1:])
        row[1:] = row[1:][sorted_row_indices]
        ICD_status = np.isin(ICD_codes, row)
        ICD_counts[i, 0] = row[0]
        ICD_counts[i, 1:] = ICD_status
        if field == "41270":
            row_status = np.isin(row[1:], ICD_codes[ICD_status])
            ICD_dates[i, 0] = row[0]
            date_row = ICD_dates_original[i][1:].astype(float)
            date_row = date_row[np.isnan(date_row) == False][sorted_row_indices]
            ICD_dates[i, 1:][ICD_status] = date_row[row_status]
    for code in HF_ICD_codes:
        is_code = subset == code
        HF_ICD_counts["counts_" + code] = np.any(is_code.to_numpy(), axis = 1)
        if field == "41270":
            case_rows = np.any(is_code.to_numpy(), axis = 1)
            HF_ICD_dates["date_" + code] = np.nan
            HF_ICD_dates.loc[case_rows, "date_" + code] = ICD_dates_original[is_code]  
    ICD_colnames = ["eid"] + ["counts_" + code for code in ICD_codes]
    ICD_counts = pd.DataFrame(ICD_counts)
    ICD_counts.columns = ICD_colnames
    ICD_counts_subsets.append(ICD_counts)
    HF_ICD_counts_subsets.append(HF_ICD_counts)
    if field == "41270":
        ICD_date_colnames = ["eid"] + ["date_" + code for code in ICD_codes]
        ICD_dates = pd.DataFrame(ICD_dates)
        ICD_dates.columns = ICD_date_colnames
        ICD_counts_subsets.append(ICD_dates)
        HF_ICD_counts_subsets.append(HF_ICD_dates)

# deletes large objects from loop
del ICD_counts
del ICD_dates

# Checks that every unique person's unique ICD code has a (valid) diagnosis date
dates_missed = []
for code in ICD_codes:
    name1, name2 = "counts_" + code, "date_" + code
    dates = ICD_counts_subsets[3].loc[ICD_counts_subsets[2][name1]!= 0, name2]
    dates_missed.append(len(dates) - np.sum(dates >= 0))

if np.sum(dates_missed) != 0:
    print("exiting: some ICD codes were not assigned a date.")
    exit()

# these are created now and used later to delete several large objects
ICD_eids = ICD_counts_subsets[0]["eid"]
ICD_eids.to_csv("ICD_eids.txt", sep = "\t", header = True, index = False)

all_ICD_dates = pd.DataFrame(ICD_counts_subsets[3].to_numpy()[:, 1:])
# these are used now
set1 = ICD_counts_subsets[0].to_numpy()[:, 1:]
set2 = ICD_counts_subsets[1].to_numpy()[:, 1:]
set3 = ICD_counts_subsets[2].to_numpy()[:, 1:]
combined_set = np.logical_or(set1, set2)
combined_set = np.logical_or(combined_set, set3)
del set1
del set2
del set3
del ICD_counts_subsets

all_ICD_cols = ["counts_" + code for code in ICD_codes]
main_ICD_cols = ["counts_" + code for code in main_ICD_codes]
all_ICD_data = pd.DataFrame(combined_set)
all_ICD_data.columns = all_ICD_cols
main_ICD_data = all_ICD_data[main_ICD_cols]
main_ICD_data.to_csv("ICD_data_main.txt", sep = "\t", header = True, index = False)
del main_ICD_data
del all_ICD_data
del combined_set

names = ["40002", "40001", "41270", "41280_as_age"]
for i in range(4):
    HF_ICD_counts_subsets[i].to_csv(names[i] + "for_HF.txt", sep = "\t", header = True, index = False)

# remaining fields: "4674", "738", "41249", "41245", "41246"
emergency_data = fields.loc[:, colnames[np.isin(colfields, ["eid", "41249"])]]
emergency_fields = np.arange(2001, 2020).astype(str)
emergency_data["was_emergency"] = np.any(np.isin(emergency_data.to_numpy(), emergency_fields), axis = 1)
HF_context = emergency_data.loc[:, ["eid", "was_emergency"]]

cardiologist_attending_data = fields.loc[:, colnames[np.isin(colfields, ["eid", "41245"])]]
cardiologist_case_data = fields.loc[:, colnames[np.isin(colfields, ["eid", "41246"])]]
cardiologist_attending = np.any(cardiologist_attending_data.to_numpy() == "1070", axis = 1)
cardiology_case = np.any(cardiologist_case_data.to_numpy() == "1160", axis = 1)
HF_context["cardiac_event"] = np.logical_or(cardiologist_attending, cardiology_case)
HF_context["cardiac_emergency"] = np.logical_and(HF_context["cardiac_event"], HF_context["was_emergency"])

insurance_data = fields.loc[:, colnames[np.isin(colfields, ["eid", "4674"])]]
insurance_data.loc[insurance_data["4674-0.0"] == "-1", "4674-0.0"] = np.nan
insurance_data.loc[insurance_data["4674-0.0"] == "-3", "4674-0.0"] = np.nan
HF_context["private insurance degree"] = insurance_data.loc[:, "4674-0.0"].astype(float)
insurance_data["private insurance"] = 0
insurance_data.loc[insurance_data["4674-0.0"] == "1", "private insurance"] = 1
insurance_data.loc[insurance_data["4674-0.0"] == "2", "private insurance"] = 1
insurance_data.loc[insurance_data["4674-0.0"] == "3", "private insurance"] = 1
HF_context["private insurance"] = insurance_data.loc[:, "private insurance"]

# 129 (north) or 130 (east) coordinate, 189 townsend deprevation index 
locality_data = fields.loc[:, colnames[np.isin(colfields, ["eid", "129", "130"])]]
locality_data[locality_data["129-0.0"] == -1] = np.nan
locality_data[locality_data["130-0.0"] == -1] = np.nan
HF_context["C1"] = locality_data["129-0.0"]
HF_context["C2"] = locality_data["130-0.0"]

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

other_data = fields.loc[:, colnames[np.isin(colfields, ["eid"] + other_fields)]]
HF_context = HF_context.merge(other_data, how = "inner", on = "eid")

health_data = fields.loc[:, colnames[np.isin(colfields, ["eid", "6150"])]]
HF_context["heart attack"] = np.any(health_data.to_numpy() == "1", axis = 1)
HF_context["angina"] = np.any(health_data.to_numpy() == "2", axis = 1)
HF_context["stroke"] = np.any(health_data.to_numpy() == "3", axis = 1)
HF_context["high BP"] = np.any(health_data.to_numpy() == "4", axis = 1)

age_data = fields.loc[:, colnames[np.isin(colfields, ["eid", "21003"])]]
HF_context["age"] = age_data.loc[:, "21003-0.0"]

income_data = fields.loc[:, colnames[np.isin(colfields, ["eid", "738"])]].astype(float)
income_data.loc[income_data["738-0.0"] == -3, "738-0.0"] = np.nan
income_data.loc[income_data["738-0.0"] == -1, "738-0.0"] = np.nan
HF_context["income"] = income_data["738-0.0"]
HF_context.to_csv("HF_metadata.txt", sep = "\t", header = True, index = False)

#---------------------------------------------------------------------------------------------------------------------------------
# Here is where I carve out the set for dimensionality reduction
# Each diagnosis will be a one hot vector
# Each patient will be an array of one hot vectors
#---------------------------------------------------------------------------------------------------------------------------------
