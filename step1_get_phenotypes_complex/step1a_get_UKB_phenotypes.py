import numpy as np
import pandas as pd
import os
import pdb
from functools import reduce
from copy import deepcopy as COPY

from step1a_library import is_field
from step1a_library import average_over_reps
from step1a_library import binarize_categoricals
from step1a_library import get_custom_features
from step1a_library import impute_N

my_fields = ["eid"] + pd.read_csv("step1.0_my_phenotypes.txt", delimiter = "\t", header = 0, dtype = str)["FieldID"].to_list()
data_types = ["eid"] + pd.read_csv("step1.0_my_phenotypes.txt", delimiter = "\t", header = 0, dtype = str)["JohnGreggValueType"].to_list()
my_fields, data_types = np.array(my_fields), np.array(data_types)

path1 = "/project/UKB_moore/UKB_50978/phenotype/penn_freeze_11132019/ukb38304.csv"
path2 = "/project/UKB_moore/UKB_50978/phenotype/penn_add_08042020/ukb43023.csv"
path3 = "/project/UKB_moore/UKB_50978/phenotype/penn_add_12082020/ukb44783.csv"
path4 = "/project/UKB_moore/UKB_50978/phenotype/penn_add_12232020/ukb44910.csv"

fields1 = pd.read_csv(path1, delimiter = ',', usecols = lambda col_name: is_field(col_name, my_fields), header = 0, low_memory = False)
fields2 = pd.read_csv(path2, delimiter = ',', usecols = lambda col_name: is_field(col_name, my_fields), header = 0, low_memory = False)
fields3 = pd.read_csv(path3, delimiter = ',', usecols = lambda col_name: is_field(col_name, my_fields), header = 0, low_memory = False)
fields4 = pd.read_csv(path4, delimiter = ',', usecols = lambda col_name: is_field(col_name, my_fields), header = 0, low_memory = False)

fields = fields1.merge(fields2, how = "inner", on = "eid")
fields = fields.merge(fields3, how = "inner", on = "eid")
fields = fields.merge(fields4, how = "inner", on = "eid")

field_cols = fields.columns[fields.columns != "eid"]
field_reps = np.array([rep.split("-")[0] for rep in field_cols])

fields2, custom_features, my_fields2, data_types2 = get_custom_features(fields, field_cols, field_reps, my_fields, data_types)

field_cols2 = fields2.columns[fields2.columns != "eid"]
field_reps2 = np.array([rep.split("-")[0] for rep in field_cols2])

# The following are only used to check imputation quality:
# [20128, 26418, 26419 ,26420 ,26421 ,26422 ,26423 ,26424 ,26425 ,26426 ,26427 ,26428 ,26429 ,26430 ,26431 ,26432 ,26433 ,26434]

# why are there fewer columns than UKB listed instances??????
# 23104, 23105

fields_for_QC_only = ["93", "94", "95", "26418", "26419", "26420", "26421", "26422", "26423", "26424", "26425", "26426", "26434", "26427", "26428", "26429", "26430", "26431", "26432", "26433", "20128"]
fields_to_impute_0 = ["3859", "2976", "2966", "20487", "20488", "20490", "20495", "20497", "20498", "20521", "20523", "20524", "20526", "20527", "20528", "20529", "20530", "20531", "21024"]
fields_to_impute_higher = ["3731", "20489", "20491", "20522", "20525"]
fields_to_impute_N = np.array(fields_to_impute_higher + fields_to_impute_0)
values_of_N = [1, 4, 4, 4, 4] + len(fields_to_impute_0)*[0]
fields_to_invert = ["2976", "2966"]

reformatted_fields = custom_features
for field, data_type in zip(my_fields2[1:], data_types2[1:]):
    rep_indices = np.where(field_reps2 == field)[0]
    replicates = fields2[['eid'] + field_cols2[rep_indices].to_list()]
    if data_type in ["Integer", "Continuous"]:
        if field in fields_to_impute_N:
            if field in fields_to_invert: invert_status = True
            else: invert_status = False
            N = values_of_N[np.where(field == fields_to_impute_N)[0][0]]
            replicates, N_was_imputed = impute_N(replicates, N, invert_status)
        reformatted_field = average_over_reps(replicates, data_type)
        if field in fields_to_impute_N:
            reformatted_field[field + "-WasModeImputed"] = N_was_imputed
        reformatted_fields.append(reformatted_field)
    else:
        if field in fields_to_impute_N: 
            if field in fields_to_invert: invert_status = True
            else: invert_status = False
            N = values_of_N[np.where(field == fields_to_impute_N)[0][0]]
            replicates, N_was_imputed = impute_N(replicates, N, invert_status)
        reformatted_field = binarize_categoricals(replicates, data_type)
        if field in fields_to_impute_N:
            reformatted_field[field + "-WasModeImputed"] = N_was_imputed
        reformatted_fields.append(reformatted_field)

if not os.path.isdir("reformatted_feature_columns"):
    os.mkdir("reformatted_feature_columns") 

col_names = []
for field in reformatted_fields:
    field_cols = field.columns[field.columns != "eid"]
    missingness_col_indices = np.array([rep.split("-")[-1] in ["unknown", "secret", "missing"] for rep in field_cols])
    QC_only_col_indices = np.array([rep.split("-")[0] in fields_for_QC_only for rep in field_cols])
    indices_to_remove = np.logical_or(QC_only_col_indices, missingness_col_indices) 
    field_cols = field_cols[indices_to_remove == False]
    for col in field_cols: 
        col_names.append(col)
        field[[col]].astype(float).to_csv("reformatted_feature_columns/" + col + ".txt", sep = "\t", header = False, index = False)

col_names = pd.DataFrame(col_names)
col_names.to_csv("reformatted_feature_columns/col_names.txt", sep = "\t", header = False, index = False)

final_fields = reduce(lambda x, y: x.merge(y, how = 'inner', on = 'eid'), reformatted_fields)
final_fields.to_csv("reformatted_fields.txt", sep = "\t", header = True, index = False)
