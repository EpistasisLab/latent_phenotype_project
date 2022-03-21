import numpy as np
import pandas as pd
import os
import pdb
epistasis_getters_folder = "epistasis_getters"
chromosomes = np.arange(1, 23).astype(str).tolist() + ["X", "XY", "Y", "MT"]
rsID_info = pd.read_csv("vQTL_hits_female.txt", delimiter = "\t", dtype = str)
rsIDs = np.unique(rsID_info["rsID"].to_numpy())
if not os.path.exists(epistasis_getters_folder):
    os.mkdir(epistasis_getters_folder)
list = open("step10g_get_commands_RUN_IN_SMALL_BATCHES.txt", 'w')
getter_prefix = "epistasis_getter_intermediaries/"
if not os.path.exists("epistasis_getter_intermediaries"):
    os.mkdir("epistasis_getter_intermediaries")
for rsID in rsIDs: 

    list.write("bsub < " + getter_prefix + "step10g_get_significant_epistasis_" + rsID + ".sh\n")
    main_file = open(getter_prefix + "step10g_get_significant_epistasis_" + rsID + ".sh", 'w')
    main_file.write("#!/bin/bash\n")
    main_file.write("#BSUB -J " + getter_prefix + "step10g_get_significant_epistasis_" + rsID + "\n")
    main_file.write("#BSUB -o " + getter_prefix + "step10g_get_significant_epistasis_" + rsID + ".out\n")
    main_file.write("#BSUB -e " + getter_prefix + "step10g_get_significant_epistasis_" + rsID + ".err\n")
    main_file.write("\n")
    for i in chromosomes:
        phenotypes = np.unique(rsID_info.loc[rsID_info["rsID"] == rsID, "phenotype"])
        for p in phenotypes:
            # makes the common part of the shell script
            suffix = "/get_epistasis_" + rsID + "_chr" + i + "_Phenotype" + str(p)
            template = "#!/bin/bash\n"
            template += "#BSUB -J " + epistasis_getters_folder + suffix + "\n"
            template += "#BSUB -o " + epistasis_getters_folder + suffix + ".out\n" 
            template += "#BSUB -e " + epistasis_getters_folder + suffix + ".err\n\n" 
            template += "source activate torch_env2\n\n" 
            # opens the shell script for the ith chromosome plink file
            path = epistasis_getters_folder + suffix + ".sh"
            main_file.write("bsub < " + path + "\n")
            file = open(path, 'w')
            # writes the shell script for the ith chromosome plink file
            content = template
            phenotypes = np.unique(rsID_info.loc[rsID_info["rsID"] == rsID, "phenotype"])
            content += "python step0_significant_epistasis_getter.py --chr " 
            content += str(i) + " --rsID " + rsID + " --pheno " + str(p) + "\n"   
            file.write(content)
            file.close()
    main_file.close()
list.close()