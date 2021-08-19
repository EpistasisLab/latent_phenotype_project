import pandas as pd
import numpy as np

LD_info = pd.read_csv("UKB_samples_unrelated_pruned.ld", delim_whitespace = True, header = 0)
print("The maximum r^2 value is " + str(np.max(LD_info["R"].to_numpy()**2)))