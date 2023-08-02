#!/bin/bash
#BSUB -J python step7a_create_AE_phenotypes
#BSUB -o python step7a_create_AE_phenotypes.out
#BSUB -e python step7a_create_AE_phenotypes.err

source activate torch_env2

DIR="AE_CV_fold_networks"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

DIR="AE_fold_error_dfs"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

DIR="AE_transformers"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for j in 1000 
do
    for k in 5
    do
        for p in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45
        do
            printf '#!/bin/bash'"\n" > $DIR'/AE_transformer'$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -J '$DIR'/AE_transformer'$j'j_'$k'k_p'$p'\n' >> $DIR"/AE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -o '$DIR'/AE_transformer'$j'j_'$k'k_p'$p'.out\n' >> $DIR"/AE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -e '$DIR'/AE_transformer'$j'j_'$k'k_p'$p'.err\n' >> $DIR"/AE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -n 20'"\n" >> $DIR"/AE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -R "rusage[mem=40000MB]"'"\n" >> $DIR"/AE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf '#BSUB -M 40000MB'"\n" >> $DIR"/AE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf 'source activate torch_env2'"\n\n" >> $DIR"/AE_transformer"$j'j_'$k'k_p'$p'.sh'
            printf 'python step7a_create_AE_phenotypes.py --rs 0 --d1 '$j' --dop '$p >> $DIR"/AE_transformer"$j'j_'$k'k_p'$p'.sh'
        done
    done
done

for j in 1000
do
    for k in 5
    do
        for p in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45
        do
            bsub < $DIR"/AE_transformer"$j'j_'$k'k_p'$p'.sh'
        done
    done
done