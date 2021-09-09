import pandas as pd
import numpy as np
from copy import deepcopy as COPY
import pdb

fam_file_name = "../step4_remove_relatives/UKB_samples_unrelated.fam"
eids = pd.read_csv(fam_file_name, delim_whitespace = True, header = None, usecols = [0], dtype = str)
eids.to_csv("eids_imputed.tab", sep = "\t", header = False, index = False)
eids[[0, 0]].to_csv("eids_real.tab", sep = "\t", header = False, index = False)
