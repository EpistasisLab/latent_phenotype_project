#!/bin/bash
#BSUB -J step7i_create_best_phenotypes_normal_AE_0.3dop
#BSUB -o step7i_create_best_phenotypes_normal_AE_0.3dop.out
#BSUB -e step7i_create_best_phenotypes_normal_AE_0.3dop.err
#BSUB -n 20
#BSUB -R "rusage[mem=40000MB]"
#BSUB -M 40000MB
source activate torch_env2

DIR="normalAE_final_model_fold_networks"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

DIR="normalAE_final_model_error_df"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

python step7i_create_best_phenotypes_normal_AE_0.3dop.py --rs 0 --ws 5 --d1 1000 --dop 0.3