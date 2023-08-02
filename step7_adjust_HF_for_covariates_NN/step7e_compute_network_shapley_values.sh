#!/bin/bash
#BSUB -J step7e_compute_network_shapley_values
#BSUB -o step7e_compute_network_shapley_values.out
#BSUB -e step7e_compute_network_shapley_values.err
source activate torch_env2

DIR="final_model_shapley_values"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

DIR="final_model_shapley_value_computers"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
do
    printf '#!/bin/bash'"\n" > $DIR'/shapley_value_computer'$i'.sh'
    printf '#BSUB -J '$DIR'/shapley_value_computer'$i'\n' >> $DIR'/shapley_value_computer'$i'.sh'
    printf '#BSUB -o '$DIR'/shapley_value_computer'$i'.out\n' >> $DIR'/shapley_value_computer'$i'.sh'
    printf '#BSUB -e '$DIR'/shapley_value_computer'$i'.err\n' >> $DIR'/shapley_value_computer'$i'.sh'
    printf '#BSUB -n 20'"\n" >> $DIR'/shapley_value_computer'$i'.sh'
    printf '#BSUB -R "rusage[mem=10000MB]"'"\n" >> $DIR'/shapley_value_computer'$i'.sh'
    printf '#BSUB -M 10000MB'"\n" >> $DIR'/shapley_value_computer'$i'.sh'
    printf 'source activate torch_env2'"\n\n" >> $DIR'/shapley_value_computer'$i'.sh'
    printf 'python step7e_compute_network_shapley_values.py --index '$i >> $DIR'/shapley_value_computer'$i'.sh'
done

for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
do
    bsub < $DIR'/shapley_value_computer'$i'.sh'
done