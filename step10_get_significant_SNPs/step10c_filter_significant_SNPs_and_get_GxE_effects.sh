#!/bin/bash
#BSUB -J step10c_filter_significant_SNPs_and_get_vQTLxE_effects
#BSUB -o step10c_filter_significant_SNPs_and_get_vQTLxE_effects.out
#BSUB -e step10c_filter_significant_SNPs_and_get_vQTLxE_effects.err

source activate torch_env2

DIR="hit_getters"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for i in {1..26}
do
    printf '#!/bin/bash'"\n" > "hit_getters/get_chr"$i"_hits.sh"
    printf '#BSUB -J step10c_filter_significant_SNPs_and_get_vQTLxE_effects'$i"\n" >> "hit_getters/get_chr"$i"_hits.sh"
    printf '#BSUB -o hit_getters/get_chr'$i"_hits.out\n" >> "hit_getters/get_chr"$i"_hits.sh"
    printf '#BSUB -e hit_getters/get_chr'$i"_hits.err\n" >> "hit_getters/get_chr"$i"_hits.sh"
    if [ "$i" == "8" ]; then
        printf '#BSUB -R "rusage[mem=50000MB]"'"\n" >> "hit_getters/get_chr"$i"_hits.sh"
        printf '#BSUB -M 50000MB'"\n" >> "hit_getters/get_chr"$i"_hits.sh"
    elif [ "$i" == "6" ]; then
        printf '#BSUB -R "rusage[mem=50000MB]"'"\n" >> "hit_getters/get_chr"$i"_hits.sh"
        printf '#BSUB -M 50000MB'"\n" >> "hit_getters/get_chr"$i"_hits.sh"
    else
        printf '#BSUB -R "rusage[mem=30000MB]"'"\n" >> "hit_getters/get_chr"$i"_hits.sh"
        printf '#BSUB -M 30000MB'"\n" >> "hit_getters/get_chr"$i"_hits.sh"
    fi
    printf 'source activate torch_env2'"\n\n" >> "hit_getters/get_chr"$i"_hits.sh"
    
    for NAME in "_alcohol" "_smoking" "_gender" "_exercise" "_age" "_BMI" "_pol"
    do
        printf 'python step0_filter_significant_SNPs_and_get_vQTLxE_effects.py --chr '$i" --name "$NAME"\n">> "hit_getters/get_chr"$i"_hits.sh"
    done
done

for i in {1..26}
do 
    bsub < hit_getters/get_chr"$i"_hits.sh
done
