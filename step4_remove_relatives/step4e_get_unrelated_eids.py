import numpy as np
import pandas as pd
import pdb
from copy import deepcopy as COPY
from tqdm import tqdm
from itertools import combinations
from matplotlib import pyplot as plt

def get_eids_to_discard(related_eids):
    
    print("step 1 complete")
    related_eids_column0 = related_eids[0].to_numpy()
    related_eids_column1 = related_eids[1].to_numpy()
    all_unique_eids = np.unique(related_eids[0])

    print("step 2 complete")  
    eid_relatives = {}
    for eid in all_unique_eids: eid_relatives[eid] = []
    for i in range(len(related_eids)): eid_relatives[related_eids_column0[i]].append(related_eids_column1[i])
    kin_counts = np.array([len(eid_relatives[eid]) for eid in all_unique_eids])

    print("step 3 complete")
    sorted_kin_count_indices = np.argsort(kin_counts)
    sorted_kin_counts = kin_counts[sorted_kin_count_indices]
    sorted_unique_eids = all_unique_eids[sorted_kin_count_indices]

    print("step 4 complete")
    num_chunks = 1000
    chunk_size = int(len(all_unique_eids)/num_chunks)
    eid_relatives2 = {}

    print("step 5 complete")
    for eid in all_unique_eids: eid_relatives2[eid] = []

    print("step 6 starting")
    for i in tqdm(range(num_chunks)):

        if i < (num_chunks - 1):
            unique_eids = sorted_unique_eids[i*chunk_size:(i+1)*chunk_size]
        else:
            unique_eids = sorted_unique_eids[i*chunk_size:]

        unique_eid_index_sets = (related_eids_column0 == unique_eids.reshape(-1,1)).astype(np.bool_)
    
        for eid, index_set in zip(unique_eids, unique_eid_index_sets):
            if len(eid_relatives2[eid]) == 0:
                relatives_to_remove = np.unique(related_eids[1][index_set])
                for relative in relatives_to_remove:
                    eid_relatives2[relative].append(eid)

    eids_to_discard = np.array([eid for eid in sorted_unique_eids if len(eid_relatives2[eid]) > 0])
    return(eids_to_discard)

print("step 1 started")
within_subjobs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
between_subjobs = list(combinations(within_subjobs, 2))
within_subjob_paths = ["kinship_within" + job + ".kin0" for job in within_subjobs]
between_subjob_paths = ["kinship_between" + job[0] + job[1] + ".kin0" for job in between_subjobs]
subjob_paths = within_subjob_paths + between_subjob_paths
related_eids1 = pd.concat([pd.read_csv(path, delim_whitespace = True, usecols = ["FID1", "FID2"], header = 0) for path in subjob_paths])
related_eids2 = COPY(related_eids1[["FID2", "FID1"]])
related_eids1.columns = [0, 1]
related_eids2.columns = [0, 1]

all_related_eids = pd.concat([related_eids1, related_eids2])
status = pd.read_csv("../step1_get_phenotypes_simple/y.txt", delimiter = "\t", header = 0, usecols = ["eid", "any_HF"])
status.columns = [0, 2] 
all_related_eids = all_related_eids.merge(status, on = 0, how = "inner")
status.columns = [1, 3] 
all_related_eids = all_related_eids.merge(status, on = 1, how = "inner")

# case eids are only excluded for relatedness to other case eids
# control eids are automatically excluded if they are related to any case eid
# control eids are also excluded for relatedness to each other
case_eid_indices = np.all(all_related_eids[[2, 3]], axis = 1)
case_eids = all_related_eids.loc[case_eid_indices, [0,1]]
control_eids = all_related_eids.loc[all_related_eids[2] == False, [0,1,3]]
case_related_control_eids = control_eids.loc[control_eids[3], [0,1]]
case_unrelated_control_eids = control_eids.loc[control_eids[3] == False, [0,1]]

case_eids_to_discard = get_eids_to_discard(case_eids)
control_eids_to_discard1 = get_eids_to_discard(case_unrelated_control_eids)
control_eids_to_discard2 = np.unique(case_related_control_eids[0].to_numpy())
eids_to_discard = np.concatenate([case_eids_to_discard, control_eids_to_discard1, control_eids_to_discard2])

all_filtered_eids = pd.read_csv("UKB_samples_filtered.fam", delimiter = " ", header = None)[0].to_numpy()
unrelated_eids = np.setdiff1d(all_filtered_eids, eids_to_discard)
if len(unrelated_eids)%10 != 0:
    remainder = len(unrelated_eids)%10
    most_control_eids = case_unrelated_control_eids[0].to_numpy()
    most_kept_control_eids = np.intersect1d(most_control_eids, unrelated_eids)
    remainder_eids = most_kept_control_eids[:remainder]
    unrelated_eids = np.setdiff1d(unrelated_eids, remainder_eids)
output = pd.DataFrame(np.array([unrelated_eids, unrelated_eids]).transpose())
output.to_csv("unrelated_eids.tab", sep = "\t", header = False, index = False)