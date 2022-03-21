#!/bin/bash
#BSUB -J step9d_get_LD_scores
#BSUB -o step9d_get_LD_scores.out
#BSUB -e step9d_get_LD_scores.err

DIR="LD_getters"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for i in {1..22} "MT" "X" "XY" "Y"
do
    printf '#!/bin/bash'"\n" > $DIR'/get_chr'$i'_LD.sh'
    printf '#BSUB -J '$DIR'/get_chr'$i'_LD'"\n" >> $DIR"/get_chr"$i"_LD.sh"
    printf '#BSUB -o '$DIR'/get_chr'$i'_LD.out'"\n" >> $DIR"/get_chr"$i"_LD.sh"
    printf '#BSUB -e '$DIR'/get_chr'$i'_LD.err'"\n" >> $DIR"/get_chr"$i"_LD.sh"
    printf '#BSUB -R "rusage[mem=10000MB]"'"\n" >> $DIR"/get_chr"$i"_LD.sh"
    printf '#BSUB -M 10000MB'"\n" >> $DIR"/get_chr"$i"_LD.sh"
    printf 'source activate torch_env2'"\n\n" >> $DIR"/get_chr"$i"_LD.sh"
    printf 'python step0_get_LD_scores.py --chr '$i >> $DIR"/get_chr"$i"_LD.sh"
done

for i in {1..22} "MT" "X" "XY" "Y"
do
    bsub < $DIR'/get_chr'$i'_LD.sh'
done