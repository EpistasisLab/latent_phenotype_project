import numpy as np
import pandas as pd
import os
import pdb
from bed_reader import open_bed
from tqdm import tqdm
from functools import reduce\

QTL_getters_folder = "QTL_getters"
chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
if not os.path.exists(QTL_getters_folder):
    os.mkdir(QTL_getters_folder)
main_file = open("step9c_get_significant_QTLs.sh", 'w')
main_file.write("#!/bin/bash\n")
main_file.write("#BSUB -J step9c_get_significant_QTLs\n")
main_file.write("#BSUB -o step9c_get_significant_QTLs.out\n")
main_file.write("#BSUB -e step9c_get_significant_QTLs.err\n")
    

end_part = ""
for i in chromosomes:
    main_file.write("\n")

    for j in range(17): 

        # makes the common part of the shell script
        template = "#!/bin/bash\n"
        template += "#BSUB -J " + QTL_getters_folder + "/get_QTLs_chr" + i + "_P" + str(j) + "\n"
        template += "#BSUB -o " + QTL_getters_folder + "/get_QTLs_chr" + i + "_P" + str(j) + ".out\n" 
        template += "#BSUB -e " + QTL_getters_folder + "/get_QTLs_chr" + i + "_P" + str(j) + ".err\n\n" 

        template += "source activate torch_env2\n\n" 

        # opens the shell script for the ith chromosome plink file
        path = QTL_getters_folder + "/get_QTLs_chr" + i + "_P" + str(j) + ".sh"
        if j == 16:
            end_part += "bsub < " + path + "\n"
        else:
            main_file.write("bsub < " + path + "\n")
        file = open(path, 'w')

        # writes the shell script for the ith chromosome plink file
        content = template
        content += "python step0_significant_QTL_getter.py --chr " + str(i) + " --pheno " + str(j) + "\n"   
        file.write(content)
        file.close()

main_file.write("\n\n" + end_part)
main_file.close()