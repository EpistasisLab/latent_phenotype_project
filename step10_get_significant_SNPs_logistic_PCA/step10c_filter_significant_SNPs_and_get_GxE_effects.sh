#!/bin/bash
#BSUB -J step10c_filter_significant_SNPs_and_get_GxE_effects
#BSUB -o step10c_filter_significant_SNPs_and_get_GxE_effects.out
#BSUB -e step10c_filter_significant_SNPs_and_get_GxE_effects.err

source activate torch_env2

DIR="hit_getters"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for i in {1..26}
do
    for NAME in "_alcohol" "_smoking" "_gender" "_exercise"
    do

    printf '#!/bin/bash'"\n" > "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf '#BSUB -J step10c_filter_significant_SNPs_and_get_GxE_effects'$i"_env"$NAME"\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf '#BSUB -o hit_getters/get_chr'$i"_hits"$NAME".out\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf '#BSUB -e hit_getters/get_chr'$i"_hits"$NAME".err\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf '#BSUB -R "rusage[mem=100000MB]"'"\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf '#BSUB -M 100000MB'"\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    printf 'source activate torch_env2'"\n\n" >> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    
    printf 'python step0_filter_significant_SNPs_and_get_GxE_effects.py --chr '$i" --name "$NAME"\n">> "hit_getters/get_chr"$i"_hits"$NAME".sh"
    done
done

for i in {1..26}
do
    for NAME in "_alcohol" "_smoking" "_gender" "_exercise"
    do 
    bsub < hit_getters/get_chr"$i"_hits"$NAME".sh
    done
done

# for i in {1..100}
# do 
#     bkill $((14597843+$i))
# done


