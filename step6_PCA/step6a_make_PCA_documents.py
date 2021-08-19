import pandas as pd
import pdb

ind_file = pd.read_csv("UKB_samples_unrelated_pruned.fam", delim_whitespace = True, header = None, usecols = [0,4,5])
ind_file.loc[:, 5] = "all_eids"
ind_file.to_csv("UKB_samples_unrelated_pruned.ind", sep = " ", header = False, index = False)

training_pop_file = open("UKB_samples_unrelated_pruned.txt", "w")
training_pop_file.write("all_eids")
training_pop_file.close()