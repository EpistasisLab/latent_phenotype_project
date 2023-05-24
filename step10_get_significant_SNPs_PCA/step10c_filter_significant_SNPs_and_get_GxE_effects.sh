#!/bin/bash
#BSUB -J step10c_filter_significant_SNPs_and_get_GxE_effects
#BSUB -o step10c_filter_significant_SNPs_and_get_GxE_effects.out
#BSUB -e step10c_filter_significant_SNPs_and_get_GxE_effects.err

source activate torch_env2

for DIR in "hit_getters" "hits_QTL_alcohol" "hits_QTL_smoking" "hits_QTL_exercise" "hits_QTL_gender" "hits_GxE_p_vals_getters_exercise" "hits_GxE_p_vals_getters_alcohol" "hits_GxE_p_vals_getters_smoking" "hits_GxE_p_vals_getters_gender"
do
    if test ! -d $DIR; then
        mkdir "$DIR"
    fi
done

for i in {1..26}
do
    for NAME in "_smoking" "_alcohol" "_gender" "_exercise"
    do

    printf '#!/bin/bash'"\n" > "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf '#BSUB -J step10c_filter_significant_SNPs_and_get_GxE_effects'$i"_env"$NAME"\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf '#BSUB -o hit_getters/get_chr'$i"_hits"$NAME".out\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf '#BSUB -e hit_getters/get_chr'$i"_hits"$NAME".err\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf '#BSUB -R "rusage[mem=200000MB]"'"\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf '#BSUB -M 200000MB'"\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf 'source activate torch_env2'"\n\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    
    printf 'python step0_filter_significant_SNPs_and_get_GxE_effects.py --chr '$i" --name "$NAME"\n">> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    done
done

for i in {1..26}
do
    for NAME in "_smoking" "_alcohol" "_gender" "_exercise"
    do 
    bsub < hit_getters/get_chr"$i"_hits"$NAME".sh
    done
done

# for i in {1..100}
# do 
#     bkill $((14597843+$i))
# done


