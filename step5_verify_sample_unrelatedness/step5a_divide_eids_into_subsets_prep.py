import numpy as np
import pandas as pd
import pdb
import os
from matplotlib import pyplot as plt

eids_fam = pd.read_csv("../step4_remove_relatives/UKB_samples_unrelated.fam", delim_whitespace = True, usecols = [0, 1], header = None)

set_size = int(len(eids_fam)/10)
if set_size != len(eids_fam)/10:
    print("exiting: sample size not divisible by 10")
output_set1 = pd.DataFrame(eids_fam[:1*set_size])
output_set2 = pd.DataFrame(eids_fam[1*set_size:2*set_size])
output_set3 = pd.DataFrame(eids_fam[2*set_size:3*set_size])
output_set4 = pd.DataFrame(eids_fam[3*set_size:4*set_size])
output_set5 = pd.DataFrame(eids_fam[4*set_size:5*set_size])
output_set6 = pd.DataFrame(eids_fam[5*set_size:6*set_size])
output_set7 = pd.DataFrame(eids_fam[6*set_size:7*set_size])
output_set8 = pd.DataFrame(eids_fam[7*set_size:8*set_size])
output_set9 = pd.DataFrame(eids_fam[8*set_size:9*set_size])
output_set10 = pd.DataFrame(eids_fam[9*set_size:10*set_size])

folder = "data_subset_content"
if not os.path.exists(folder):
    os.mkdir(folder)
output_path1 = folder + "/subset1.tab"
output_path2 = folder + "/subset2.tab"
output_path3 = folder + "/subset3.tab"
output_path4 = folder + "/subset4.tab"
output_path5 = folder + "/subset5.tab"
output_path6 = folder + "/subset6.tab"
output_path7 = folder + "/subset7.tab"
output_path8 = folder + "/subset8.tab"
output_path9 = folder + "/subset9.tab"
output_path10 = folder + "/subset10.tab"

output_set1.to_csv(output_path1, sep = "\t", header = False, index = False)
output_set2.to_csv(output_path2, sep = "\t", header = False, index = False)
output_set3.to_csv(output_path3, sep = "\t", header = False, index = False)
output_set4.to_csv(output_path4, sep = "\t", header = False, index = False)
output_set5.to_csv(output_path5, sep = "\t", header = False, index = False)
output_set6.to_csv(output_path6, sep = "\t", header = False, index = False)
output_set7.to_csv(output_path7, sep = "\t", header = False, index = False)
output_set8.to_csv(output_path8, sep = "\t", header = False, index = False)
output_set9.to_csv(output_path9, sep = "\t", header = False, index = False)
output_set10.to_csv(output_path10, sep = "\t", header = False, index = False)
