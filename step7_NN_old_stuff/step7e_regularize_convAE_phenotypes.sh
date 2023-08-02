#!/bin/bash
#BSUB -J step7e_regularize_convAE_phenotypes
#BSUB -o step7e_regularize_convAE_phenotypes.out
#BSUB -e step7e_regularize_convAE_phenotypes.err

source activate torch_env2

DIR="dropout_CV_fold_networks"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

DIR="dropout_fold_error_dfs"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

DIR="dropout_convAE_transformers"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for j in 4000
do
    for k in 1 5
    do
        for p in 0.01 0.05 0.1 0.15 0.2 0.25
        do
            printf '#!/bin/bash'"\n" > $DIR'/convAE_transformer'$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -J '$DIR'/convAE_transformer'$j'j_'$k'k_p'$p'\n' >> $DIR"/convAE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -o '$DIR'/convAE_transformer'$j'j_'$k'k_p'$p'.out\n' >> $DIR"/convAE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -e '$DIR'/convAE_transformer'$j'j_'$k'k_p'$p'.err\n' >> $DIR"/convAE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -n 20'"\n" >> $DIR"/convAE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -R "rusage[mem=40000MB]"'"\n" >> $DIR"/convAE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -M 40000MB'"\n" >> $DIR"/convAE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf 'source activate torch_env2'"\n\n" >> $DIR"/convAE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf 'python step7e_regularize_convAE_phenotypes.py --rs 0 --ws '$k' --d1 '$j' --dop '$p >> $DIR"/convAE_transformer"$j'j_'$k'k_p'$p'.sh'
        done
    done
done

#bsub < dropout_convAE_transformers/convAE_transformer40j_1k.sh
for j in 4000
do
    for k in 1 5
    do
        for p in 0.01 0.05 0.1 0.15 0.2 0.25
        do
            bsub < $DIR"/convAE_transformer"$j'j_'$k'k_p'$p'.sh'
        done
    done
done