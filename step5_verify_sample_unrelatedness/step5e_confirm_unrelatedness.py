import numpy as np
import pandas as pd
import pdb
from copy import deepcopy as COPY
from tqdm import tqdm
from itertools import combinations
from matplotlib import pyplot as plt

within_subjobs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
between_subjobs = list(combinations(within_subjobs, 2))
within_subjob_paths = ["kinship_within" + job + ".kin0" for job in within_subjobs]
between_subjob_paths = ["kinship_between" + job[0] + job[1] + ".kin0" for job in between_subjobs]
subjob_paths = within_subjob_paths + between_subjob_paths
related_eids = pd.concat([pd.read_csv(path, delim_whitespace = True, usecols = ["FID1", "FID2"], header = 0) for path in subjob_paths])
num_related_pairs = str(len(related_eids))
print("There are " + num_related_pairs + " related pairs.")