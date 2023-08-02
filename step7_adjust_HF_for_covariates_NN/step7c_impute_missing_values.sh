#!/bin/bash
#BSUB -J step7c_impute_missing_values
#BSUB -o step7c_impute_missing_values.out
#BSUB -e step7c_impute_missing_values.err
source activate torch_env2

DIR="info_features"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

DIR="imputation_testors"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for nn in 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1
do
        printf '#!/bin/bash'"\n" > $DIR'/imputation_testor'$nn'nn.sh'
        printf '#BSUB -J '$DIR'/imputation_testor'$nn'nn_\n' >> $DIR"/imputation_testor"$nn'nn.sh'
        printf '#BSUB -o '$DIR'/imputation_testor'$nn'nn.out\n' >> $DIR"/imputation_testor"$nn'nn.sh'
        printf '#BSUB -e '$DIR'/imputation_testor'$nn'nn.err\n' >> $DIR"/imputation_testor"$nn'nn.sh'
        printf '#BSUB -n 1'"\n" >> $DIR"/imputation_testor"$nn'nn.sh'
        printf '#BSUB -R "rusage[mem=80000MB]"'"\n" >> $DIR"/imputation_testor"$nn'nn.sh'
        printf '#BSUB -M 80000MB'"\n" >> $DIR"/imputation_testor"$nn'nn.sh'
        printf 'source activate torch_env2'"\n\n" >> $DIR"/imputation_testor"$nn'nn.sh'
        printf 'python step7c_impute_missing_values.py --nn '$nn >> $DIR"/imputation_testor"$nn'nn.sh'
done

for nn in 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1
do
    bsub < $DIR"/imputation_testor"$nn'nn.sh'
done