#!/bin/bash
#BSUB -J step7b_train_w2v_models
#BSUB -o step7b_train_w2v_models.out
#BSUB -e step7b_train_w2v_models.err

source activate torch_env2

DIR="w2v_transformers"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for j in 1 2 3 4 5
do
    printf '#!/bin/bash'"\n" > $DIR'/w2v_transformer'$j'j.sh'
    printf '#BSUB -J '$DIR'/w2v_transformer'$j"j\n" >> $DIR"/w2v_transformer"$j"j.sh"
    printf '#BSUB -o '$DIR'/w2v_transformer'$j"j.out\n" >> $DIR"/w2v_transformer"$j"j.sh"
    printf '#BSUB -e '$DIR'/w2v_transformer'$j"j.err\n" >> $DIR"/w2v_transformer"$j"j.sh"
    printf '#BSUB -n 20'"\n" >> $DIR"/w2v_transformer"$j"j.sh"
    printf 'source activate torch_env2'"\n\n" >> $DIR"/w2v_transformer"$j"j.sh"
    printf 'python step7b_train_w2v_models.py --rs 0 --ws '$j >> $DIR"/w2v_transformer"$j"j.sh"
done

#bsub < w2v_transformers/w2v_transformer1j.sh
for j in 1 2 3 4 5
do
    bsub < $DIR"/w2v_transformer"$j"j.sh"
done