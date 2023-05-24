#!/bin/bash
#BSUB -J step10g_get_CV_testing_accuracy
#BSUB -o step10g_get_CV_testing_accuracy.out
#BSUB -e step10g_get_CV_testing_accuracy.err

source activate torch_env2

DIR="CV_params" 
if test ! -d $DIR; then
    mkdir "$DIR"
fi

DIR="CV_param_getters" 
if test ! -d $DIR; then
    mkdir "$DIR"
fi


for i in {0..29}
do
    printf '#!/bin/bash'"\n" > "CV_param_getters/get_CV"$i".sh"
    printf '#BSUB -J step10g_get_CV_testing_accuracy_CV'$i"\n" >> "CV_param_getters/get_CV"$i".sh"
    printf '#BSUB -o CV_param_getters/get_CV'$i".out\n" >> "CV_param_getters/get_CV"$i".sh"
    printf '#BSUB -e CV_param_getters/get_CV'$i".err\n" >> "CV_param_getters/get_CV"$i".sh"
    printf '#BSUB -R "rusage[mem=10000MB]"'"\n" >> "CV_param_getters/get_CV"$i".sh"
    printf '#BSUB -M 10000MB'"\n" >> "CV_param_getters/get_CV"$i".sh"
    printf 'source activate torch_env2'"\n\n" >> "CV_param_getters/get_CV"$i".sh"
    
    printf  'python step10g_get_CV_testing_accuracy.py --CV '$i"\n">> "CV_param_getters/get_CV"$i".sh"
done

for i in {0..29}
do
    bsub < "CV_param_getters/get_CV"$i".sh"
done
