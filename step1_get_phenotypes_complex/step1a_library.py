import numpy as np
import pandas as pd
import pdb
import statsmodels.api as sm
from functools import reduce
from copy import deepcopy as COPY
from scipy.stats import mannwhitneyu
from matplotlib import pyplot as plt

def impute_N(replicates, N, invert_status):
    col_names = replicates.columns
    field_name_indices = np.where(col_names != "eid")[0]
    replicates = replicates.to_numpy()
    nan_indices = np.isnan(replicates)
    N_was_imputed = np.any(nan_indices, axis = 1)
    if invert_status == True:
        value_indices = (nan_indices == False)
        value_indices[:, 0] = False
        special_values = np.array([-7, -17, -27, -2, -3, -13, -23, -818, -1, -11, -21, -121])
        special_value_indices = np.isin(replicates, special_values)
        value_indices[special_value_indices] = False
        replicates[value_indices] = 1/replicates[value_indices]
    replicates[nan_indices] = N
    replicates = pd.DataFrame(replicates)
    replicates.columns = col_names
    return(replicates, N_was_imputed)

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

def change_special_values(replicates, data_type):

    """
    Purpose
    -------
    Parameters
    ----------
    Returns
    -------
    """

    replicate_names = replicates.columns
    replicates = replicates.to_numpy(dtype = float)

    is_missing = np.all(np.isnan(replicates[:, replicate_names !=  "eid"]), axis = 1)

    # -7, -17, -27, and -2 are each, for all practical purposes, "None of the above". 
    replicates[replicates == -2] = -7
    replicates[replicates == -17] = -7
    replicates[replicates == -27] = -7

    # -3, -13, -23, and -818 are literally the same
    replicates[replicates == -818] = -3
    replicates[replicates == -13] = -3
    replicates[replicates == -23] = -3
    is_a_secret = np.any(replicates == -3, axis = 1)
    replicates[replicates == -3] = np.nan

    # -1, -11, -21, and -121 are literally the same
    replicates[replicates == -121] = -1
    replicates[replicates == -11] = -1
    replicates[replicates == -21] = -1
    is_unknown = np.any(replicates == -1, axis = 1)
    replicates[replicates == -1] = np.nan

    if data_type in ["Integer", "Continuous"]:
        # -10 gets its own category if it is a category, and it gets 0.5 if it is numeric 
        replicates[replicates == -10] = 0.5

    replicates = pd.DataFrame(replicates)
    replicates.columns = replicate_names
    return(replicates, is_unknown, is_a_secret, is_missing)

def average_over_reps(replicates, data_type):

    """
    Purpose
    -------
    Parameters
    ----------
    Returns
    -------
    """

    replicates, is_unknown, is_a_secret, is_missing = change_special_values(replicates, data_type)

    field_names = replicates.columns[replicates.columns != 'eid']
    field = field_names[0].split("-")[0]
    rep_values = replicates[field_names].to_numpy()
    rep_values = rep_values.reshape(-1, len(rep_values[0]))

    # at least one numeric field (26416) has some -7 values for no discernable reason
    # -7 normally means (none of the above), which is equivalent to being missing for fields over a numeric domain. 
    rep_values[rep_values == -7] = np.nan

    averages = replicates.loc[:, ["eid"]]
    averages[field + "-average"] = np.nanmean(rep_values, axis = 1)
    averages[field + "-is-a-secret"] = is_a_secret
    averages[field + "-is-unknown"] = is_unknown
    averages[field + "-is-missing"] = is_missing

    return(averages)


def binarize_categoricals(replicates, data_type, nametag = ''):

    """
    Purpose
    -------
    Parameters
    ----------
    Returns
    -------
    """

    replicates, is_unknown, is_a_secret, is_missing = change_special_values(replicates, data_type)
    field_cols = replicates.columns[replicates.columns != "eid"]
    field_values = replicates.loc[:, field_cols].to_numpy()
    field = field_cols[0].split("-")[0]
    all_missing_indices = np.all(np.isnan(field_values), axis = 1)

    unique_values = np.unique(field_values)[np.isnan(np.unique(field_values)) == False]
    field_values_bin = [np.any(field_values == value, axis = 1) for value in unique_values]
    field_values_bin = pd.DataFrame(np.array(field_values_bin).astype(int).T)
    field_values_bin.columns = [field + "-" + str(val) for val in unique_values]
    field_values_bin.loc[all_missing_indices, :] = np.nan

    if data_type == "Binary":
        field_values_bin = field_values_bin[field_values_bin.columns[:-1]]

    # -7 being True is 100% correlated to all normal values being false, making it redundant. 
    if ((field + "-" + str(-7.0)) in field_values_bin.columns):
        del field_values_bin[(field + "-" + str(-7.0))]

    field_values_bin["eid"] = replicates.loc[:, "eid"]
    field_values_bin[field + nametag + "-is-unknown"] = is_unknown
    field_values_bin[field + nametag + "-is-a-secret"] = is_a_secret
    field_values_bin[field + nametag + "-is-missing"] = is_missing

    return(field_values_bin)

'''
The following function treats these special features distinctly:
1) 20107, 20110, 20111 do not follow standard formatting
   -17, -27 = No Illnesses; -11, -21 = don't know, -13, -23 = Prefer not to answer
2) 6153 and 6177 are essentially the same and to be merged
3) 46 and 47 are left and right grip strength: choose max to avoid variance from strokes.
4) 93, 94, and 95 are manual measurements to be merged with 4080, 4079, and 102 respectively 
5) 3731 (Former alcohol drinker), while binary, is to be set to "No" for naturally missing values (current alcohol drinkers). 
6) 20118 needs to be recoded as "urban" and "not urban". 
'''

def get_custom_features(fields, field_cols, field_reps, my_fields, data_types):

    fields_to_remove = []

    pack_years_fields = ["eid", "20116-0.0", "20161-0.0"]
    fields_to_remove += ["20161"]
    pack_years = fields.loc[:, pack_years_fields]
    no_pack_year_info = np.isnan(pack_years.loc[:, "20161-0.0"])
    never_smoked = (pack_years.loc[:, "20116-0.0"] == 0)
    true_0_pack_years = np.logical_and(no_pack_year_info, never_smoked)
    pack_years.loc[true_0_pack_years, "20161-0.0"] = 0
    pack_years = pack_years[["eid", "20161-0.0"]]
    pack_years.columns = ["eid", "pack-years"]

    weekly_alcohol_fields = ["eid", "1568-0.0", "1578-0.0", "1588-0.0", "1598-0.0", "1608-0.0", "5364-0.0"]
    monthly_alcohol_fields = ["eid", "4407-0.0", "4418-0.0", "4429-0.0", "4440-0.0", "4451-0.0", "4462-0.0"]
    fields_to_remove += ["1568", "1578", "1588", "1598", "1608", "5364"]
    fields_to_remove +=["4407", "4418", "4429", "4440", "4451", "4462"]

    weekly_alcohol = fields[weekly_alcohol_fields].to_numpy()
    is_unknown = np.any(weekly_alcohol == -1, axis = 1) 
    is_a_secret = np.any(weekly_alcohol == -3, axis = 1)
    is_missing = np.all(np.isnan(weekly_alcohol[:, 1:]), axis = 1)
    weekly_alcohol[is_unknown, 1:] = np.nan
    weekly_alcohol[is_a_secret, 1:] = np.nan
    weekly_alcohol = pd.DataFrame(weekly_alcohol)
    annual_alcohol1 = weekly_alcohol.loc[:, [0]]
    annual_alcohol1["annual-consumption"] = weekly_alcohol.loc[:, [1, 2, 3, 4, 5, 6]].sum(axis = 1, min_count = 1)*52.1429

    monthly_alcohol = fields[monthly_alcohol_fields].to_numpy()
    is_unknown = np.any(monthly_alcohol == -1, axis = 1) 
    is_a_secret = np.any(monthly_alcohol == -3, axis = 1)
    is_missing = np.all(np.isnan(monthly_alcohol[:, 1:]), axis = 1)
    monthly_alcohol[is_unknown, 1:] = np.nan
    monthly_alcohol[is_a_secret, 1:] = np.nan
    monthly_alcohol = pd.DataFrame(monthly_alcohol)
    annual_alcohol2 = monthly_alcohol.loc[:, [0]]
    annual_alcohol2["annual-consumption"] = monthly_alcohol.loc[:, [1, 2, 3, 4, 5, 6]].sum(axis = 1, min_count = 1)*12

    fields_to_remove.append("1558")
    drinking_status = fields[["eid", "1558-0.0"]]
    annual_alcohol3 = drinking_status.loc[:, ["eid"]]
    annual_alcohol3["annual-consumption"] = drinking_status.loc[:, "1558-0.0"] == 6
    annual_alcohol3.loc[annual_alcohol3["annual-consumption"] == 0, "annual-consumption"] = np.nan
    annual_alcohol3.loc[annual_alcohol3["annual-consumption"] == 1, "annual-consumption"] = 0

    annual_alcohol = COPY(annual_alcohol1.loc[:, :])
    nan_indices = np.isnan(annual_alcohol["annual-consumption"].to_numpy())
    annual_alcohol.loc[nan_indices, "annual-consumption"] = annual_alcohol2.loc[nan_indices, "annual-consumption"]
    nan_indices = np.isnan(annual_alcohol["annual-consumption"].to_numpy())
    annual_alcohol.loc[nan_indices, "annual-consumption"] = annual_alcohol3.loc[nan_indices, "annual-consumption"]
    annual_alcohol.columns = ["eid", "annual-consumption"]

    fields_to_remove += ["6153", "6177"]
    field6153 = fields[["eid"] + field_cols[field_reps == "6153"].to_list()]
    field6153 = binarize_categoricals(field6153, "Categorical multiple")
    field6177 = fields[["eid"] + field_cols[field_reps == "6177"].to_list()]
    field6177 = binarize_categoricals(field6177, "Categorical multiple")

    nan_indices_6177 = np.isnan(field6177["6177-1.0"])
    medications = COPY(field6177.loc[:, :])
    medications["6177-4.0"] = 0.0
    medications["6177-5.0"] = 0.0
    ordered_cols = ['6177-1.0', '6177-2.0', '6177-3.0', '6177-4.0', '6177-5.0', 'eid', '6177-is-unknown', '6177-is-a-secret', '6177-is-missing']
    medications = medications.loc[:, ordered_cols] 
    medications[nan_indices_6177] = field6153.loc[nan_indices_6177].to_numpy() 

    fields_to_remove.append("20118")
    urban_values = [1.0, 5.0, 11.0, 12.0]
    field20118 = fields[["eid"] + field_cols[field_reps == "20118"].to_list()]
    field20118 = binarize_categoricals(field20118, "Categorical single")
    urban_or_not = COPY(field20118.loc[:, ["eid"]])
    urban_or_not["is_urban"] = field20118.loc[:, ["20118-" + str(i) for i in urban_values]].sum(axis = 1)

    fields_to_remove += ["93", "94", "95"]
    fields_to_remove += ["4080", "4079", "102"]

    field_set1 = [fields[["eid"] + field_cols[field_reps == i].to_list()] for i in ["93", "94", "95"]]
    field_set1 = [average_over_reps(set, "Continuous") for set in field_set1]
    field_set2 = [fields[["eid"] + field_cols[field_reps == i].to_list()] for i in ["4080", "4079", "102"]]
    field_set2 = [average_over_reps(set, "Continuous") for set in field_set2]

    field_set = []
    for field1, field2, name1, name2 in zip(field_set1, field_set2, ["93", "94", "95"], ["4080", "4079", "102"]):
        nan_indices = field2[field2.columns[1]].isna()
        field2.loc[nan_indices, field2.columns[1]] = field1.loc[nan_indices, field1.columns[1]]
        field2.loc[nan_indices, field2.columns[2]] = field1.loc[nan_indices, field1.columns[2]]
        field2.loc[nan_indices, field2.columns[3]] = field1.loc[nan_indices, field1.columns[3]]
        field2.loc[nan_indices, field2.columns[4]] = field1.loc[nan_indices, field1.columns[4]]
        remaining_nan_indices = field2[field2.columns[1]].isna()
        manual_measurements = np.logical_and(nan_indices, remaining_nan_indices == False)
        field2["substituted-" + name1 + "WITH" + name2] = manual_measurements
        field_set.append(field2)

    bp_and_pulse = reduce(lambda x, y: x.merge(y, how = 'inner', on = 'eid'), field_set)

    # note: no instances of Breast cancer (5) exist for this field despite the website saying so. There is also no (7).
    # Group 1 : Heart disease (1), Stroke (2), High blood pressure (8), Chronic bronchitis/emphysema (6), Alzheimer's disease/dementia (10), Diabetes (9).
    # Group 2 : Parkinson's disease (11), Severe Depression (12), Lung cancer (3), Bowel cancer (4), Prostate cancer (13).
    # see https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20107 for more information 

    fields_to_remove += ["20107", "20110", "20111"]
    # replacing -13 and -23 with -11 and -21 respectively
    group1_unique_vals = [-17, -13, -11, 1, 2, 6, 8, 9, 10]
    group2_unique_vals = [-27, -23, -21, 3, 4, 5, 11, 12, 13]

    field_subgroups = []
    for field in ["20107", "20110", "20111"]:
       
        col_names = fields[field_cols[field_reps == field]].columns
        df_group1 = fields[field_cols[field_reps == field]].to_numpy()
        df_group2 = COPY(df_group1)

        df_group1[np.isin(df_group1, group2_unique_vals)] = np.nan
        df_group2[np.isin(df_group2, group1_unique_vals)] = np.nan

        df_group1, df_group2 = pd.DataFrame(df_group1), pd.DataFrame(df_group2)
        df_group1.columns = col_names
        df_group2.columns = col_names
        df_group1["eid"], df_group2["eid"] = fields.loc[:, "eid"], fields.loc[:, "eid"]
        
        field_subgroups += [df_group1, df_group2]

    nametags = ["-group1", "-group2", "-group1", "-group2", "-group1", "-group2"]
    binarized_dataframes = [binarize_categoricals(subgroup, "Categorical multiple", nametag = tag) for subgroup, tag in zip(field_subgroups, nametags)]
    family_history = reduce(lambda x, y: x.merge(y, how = 'inner', on = 'eid'), binarized_dataframes)

    fields_to_remove += ["46", "47"]
    fields_46_47 = fields[field_cols[np.isin(field_reps, ["46", "47"])]].to_numpy()
    max_grip_strength = fields.loc[:, ["eid"]]
    max_grip_strength["max-46-47"] = np.max(fields_46_47, axis = 1)

    custom_features = [pack_years, medications, urban_or_not, annual_alcohol, bp_and_pulse, family_history, max_grip_strength]

    columns_to_remove = field_cols[np.isin(field_reps, fields_to_remove)] 
    columns_to_keep = np.setdiff1d(fields.columns, columns_to_remove)
    fields = fields.loc[:, columns_to_keep]

    fields_to_keep = np.isin(my_fields, np.array(fields_to_remove)) == False
    my_fields, data_types = my_fields[fields_to_keep], data_types[fields_to_keep]
 
    return(fields, custom_features, my_fields, data_types)
