#!/bin/bash
#BSUB -J step10c_filter_significant_SNPs_and_get_GxE_effects
#BSUB -o step10c_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.out
#BSUB -e step10c_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.err

source activate torch_env2

for DIR in "hits_GxE_p_vals" "hits_GxE_p_vals_getters_exercise" "hits_GxE_p_vals_getters_alcohol" "hits_GxE_p_vals_getters_smoking" "hits_GxE_p_vals_getters_gender"
do
    if test ! -d $DIR; then
        mkdir "$DIR"
    fi
done

for i in {1..26}
do
    for NAME in "_smoking" "_alcohol" "_gender" "_exercise"
    do
        python step0_filter_significant_SNPs_and_get_GxE_effects_LAST_PART_ONLY.py --chr "$i" --name "$NAME"
    done
done

