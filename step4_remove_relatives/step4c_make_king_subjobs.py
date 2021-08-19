import numpy as np
import pandas as pd
import os
from itertools import combinations

within_subjobs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
between_subjobs = combinations(within_subjobs, 2)
subjob_size = len(pd.read_csv("data_subset_content/subset1.tab", delimiter = "\t", header = None))

folder = "king_subjobs"
if not os.path.exists(folder):
    os.mkdir(folder)

for job in within_subjobs:

    file = open("king_subjobs/get_kinships_" + job + ".sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#BSUB -J king\n")
    file.write("#BSUB -o king_subjobs/get_kinships_" + job + ".out\n")
    file.write("#BSUB -e king_subjobs/get_kinships_" + job + ".error\n")
    file.write('#BSUB -R "rusage[mem=40000MB]"\n')
    file.write("#BSUB -M 40000MB\n")
    file.write("#BSUB -n 5\n")
    file.write("module load KING/2.2.4\n\n")

    file.write("king -b data_subset_content/UKB_samples_filtered_subset" + job)
    file.write(".bed --kinship --degree 3 --cpus 5 --prefix kinship_within" + job)
    file.close()

for subset1, subset2 in between_subjobs:
    
    file = open("king_subjobs/get_kinships_" + subset1 + subset2 + ".sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#BSUB -J king\n")
    file.write("#BSUB -o king_subjobs/get_kinships_" + subset1 + subset2 + ".out\n")
    file.write("#BSUB -e king_subjobs/get_kinships_" + subset1 + subset2 + ".error\n")
    file.write('#BSUB -R "rusage[mem=40000MB]"\n')
    file.write("#BSUB -M 40000MB\n")
    file.write("#BSUB -n 5\n")
    file.write("module load KING/2.2.4\n\n")

    file.write("king -b data_subset_content/UKB_samples_filtered_subset" + subset1 + ".bed,")
    file.write("data_subset_content/UKB_samples_filtered_subset" + subset2 + ".bed")
    file.write(" --kinship --degree 3 --proj " + str(subjob_size))
    file.write(" --cpus 5 --prefix kinship_between" + subset1 + subset2)