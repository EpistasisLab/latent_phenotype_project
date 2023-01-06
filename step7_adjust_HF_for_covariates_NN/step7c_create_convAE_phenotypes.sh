#!/bin/bash
#BSUB -J step7c_create_convAE_phenotypes
#BSUB -o step7c_create_convAE_phenotypes.out
#BSUB -e step7c_create_convAE_phenotypes.err

source activate torch_env2

DIR="CV_fold_networks"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

DIR="fold_error_dfs"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

DIR="convAE_transformers"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for j in 1000 2000 3000 4000 4500 5000 5500
do
    for k in 1 2 3 4 5
    do
        printf '#!/bin/bash'"\n" > $DIR'/convAE_transformer'$j'j_'$k'k.sh'
        printf '#BSUB -J '$DIR'/convAE_transformer'$j'j_'$k'k'"\n" >> $DIR"/convAE_transformer"$j'j_'$k"k.sh"
        printf '#BSUB -o '$DIR'/convAE_transformer'$j'j_'$k'k.out'"\n" >> $DIR"/convAE_transformer"$j'j_'$k"k.sh"
        printf '#BSUB -e '$DIR'/convAE_transformer'$j'j_'$k'k.err'"\n" >> $DIR"/convAE_transformer"$j'j_'$k"k.sh"
        printf '#BSUB -n 20'"\n" >> $DIR"/convAE_transformer"$j'j_'$k"k.sh"
        printf '#BSUB -R "rusage[mem=40000MB]"'"\n" >> $DIR"/convAE_transformer"$j'j_'$k"k.sh"
        printf '#BSUB -M 40000MB'"\n" >> $DIR"/convAE_transformer"$j'j_'$k"k.sh"
        printf 'source activate torch_env2'"\n\n" >> $DIR"/convAE_transformer"$j'j_'$k"k.sh"
        printf 'python step7c_create_convAE_phenotypes.py --rs 0 --ws '$k' --d1 '$j >> $DIR"/convAE_transformer"$j'j_'$k"k.sh"
    done
done

#bsub < convAE_transformers/convAE_transformer40j_1k.sh
for j in 1000 2000 3000 4000 4500 5000 5500
do
    for k in 1 2 3 4 5
    do
        bsub < $DIR"/convAE_transformer"$j'j_'$k"k.sh"
    done
done