import numpy as np
import pandas as pd
import os
import pdb

vQTL_getters_folder = "vQTL_getters"
chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]
chromosomes += ["14", "15", "16", "17", "18", "19", "20", "21", "22", "MT", "X", "XY", "Y"]
if not os.path.exists(vQTL_getters_folder):
    os.mkdir(vQTL_getters_folder)
main_file = open("step9d_get_significant_vQTLs.sh", 'w')
main_file.write("#!/bin/bash\n")
main_file.write("#BSUB -J step9d_get_significant_vQTLs\n")
main_file.write("#BSUB -o step9d_get_significant_vQTLs.out\n")
main_file.write("#BSUB -e step9d_get_significant_vQTLs.err\n")


for i in chromosomes:
    main_file.write("\n")
    for j in range(16): 

        # makes the common part of the shell script
        template = "#!/bin/bash\n"
        template += "#BSUB -J " + vQTL_getters_folder + "/get_vQTLs_chr" + i + "_P" + str(j) + "\n"
        template += "#BSUB -o " + vQTL_getters_folder + "/get_vQTLs_chr" + i + "_P" + str(j) + ".out\n" 
        template += "#BSUB -e " + vQTL_getters_folder + "/get_vQTLs_chr" + i + "_P" + str(j) + ".err\n\n" 

        template += "source activate torch_env2\n\n" 

        # opens the shell script for the ith chromosome plink file
        path = vQTL_getters_folder + "/get_vQTLs_chr" + i + "_P" + str(j) + ".sh"
        main_file.write("bsub < " + path + "\n")
        file = open(path, 'w')

        # writes the shell script for the ith chromosome plink file
        content = template
        content += "python step0_significant_vQTL_getter.py --chr " + str(i) + " --pheno " + str(j) + "\n"   
        file.write(content)
        file.close()