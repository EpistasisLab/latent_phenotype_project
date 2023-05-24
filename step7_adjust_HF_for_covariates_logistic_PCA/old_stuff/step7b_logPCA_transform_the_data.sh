#!/bin/bash
#BSUB -J step7b_logPCA_transform_the_data
#BSUB -o step7b_logPCA_transform_the_data.out
#BSUB -e step7b_logPCA_transform_the_data.err

DIR="logPCA_transformers"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for i in {1..20}
do
    printf '#!/bin/bash'"\n" > $DIR'/logPCA_transformer'$i'.sh'
    printf '#BSUB -J '$DIR'/logPCA_transformer'$i'_LD'"\n" >> $DIR"/logPCA_transformer"$i".sh"
    printf '#BSUB -o '$DIR'/logPCA_transformer'$i'_LD.out'"\n" >> $DIR"/logPCA_transformer"$i".sh"
    printf '#BSUB -e '$DIR'/logPCA_transformer'$i'_LD.err'"\n" >> $DIR"/logPCA_transformer"$i".sh"
    printf '#BSUB -R "rusage[mem=10000MB]"'"\n" >> $DIR"/logPCA_transformer"$i".sh"
    printf '#BSUB -M 10000MB'"\n" >> $DIR"/logPCA_transformer"$i".sh"
    printf 'source activate torch_env2'"\n\n" >> $DIR"/logPCA_transformer"$i".sh"
    printf 'python step7b_logPCA_transform_the_data.py --k '$i >> $DIR"/logPCA_transformer"$i".sh"
done

for i in {1..20}
do
    bsub < $DIR'/logPCA_transformer'$i'.sh'
done
