#!/bin/bash
#BSUB -J step7b_create_best_phenotypes_normal_AE_0.3dop
#BSUB -o step7b_create_best_phenotypes_normal_AE_0.3dop.out
#BSUB -e step7b_create_best_phenotypes_normal_AE_0.3dop.err
#BSUB -n 20
#BSUB -R "rusage[mem=40000MB]"
#BSUB -M 40000MB
source activate torch_env2

DIR="AE_final_model_network"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

DIR="AE_final_model_error_df"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

python step7b_create_best_phenotypes_normal_AE_0.3dop.py --rs 0 --d1 1000 --dop 0.3