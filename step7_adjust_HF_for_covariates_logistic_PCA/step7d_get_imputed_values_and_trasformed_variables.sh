#!/bin/bash
#BSUB -J step7d_get_imputed_values_and_trasformed_variables
#BSUB -o step7d_get_imputed_values_and_trasformed_variables.out
#BSUB -e step7d_get_imputed_values_and_trasformed_variables.err
source activate torch_env2

python step7d_get_imputed_values_and_trasformed_variables.py