#!/bin/bash
#BSUB -J step10d_get_significant_GxE_p_values
#BSUB -o step10d_get_significant_GxE_p_values.out
#BSUB -e step10d_get_significant_GxE_p_values.err

source activate torch_env2

FILES1="/home/greggj/pleiotropy_and_GxE/step10_get_significant_SNPs_PCA/hits_GxE_p_vals_getters_alcohol/*.sh"
FILES2="/home/greggj/pleiotropy_and_GxE/step10_get_significant_SNPs_PCA/hits_GxE_p_vals_getters_exercise/*.sh"
FILES3="/home/greggj/pleiotropy_and_GxE/step10_get_significant_SNPs_PCA/hits_GxE_p_vals_getters_gender/*.sh"
FILES4="/home/greggj/pleiotropy_and_GxE/step10_get_significant_SNPs_PCA/hits_GxE_p_vals_getters_smoking/*.sh"


for f in $FILES1
do 
    bsub < $f
done

for f in $FILES2
do 
    bsub < $f
done

for f in $FILES3
do 
    bsub < $f
done

for f in $FILES4
do 
    bsub < $f
done