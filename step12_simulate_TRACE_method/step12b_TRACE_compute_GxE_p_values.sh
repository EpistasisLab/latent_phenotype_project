#!/bin/bash
#BSUB -J step12b_TRACE_compute_GxE_p_values
#BSUB -o step12b_TRACE_compute_GxE_p_values.out
#BSUB -e step12b_TRACE_compute_GxE_p_values.err

rm -f step0_GxE_p_value_getters/*.sh

for j in $(seq 0 4);
do
    for i in $(seq 0 1999); 
    do 
        touch "step0_GxE_p_value_getters/${j}_${i}.sh"
        printf '#!/bin/bash'"\n" > "step0_GxE_p_value_getters/${j}_${i}.sh"
        printf "#BSUB -J step0_GxE_p_value_getters/step12b\n" >> "step0_GxE_p_value_getters/${j}_${i}.sh"
        printf "#BSUB -o step0_GxE_p_value_getters/step12b.out\n" >> "step0_GxE_p_value_getters/${j}_${i}.sh"
        printf "#BSUB -e step0_GxE_p_value_getters/step12b.err\n" >> "step0_GxE_p_value_getters/${j}_${i}.sh"
        printf "source activate torch_env2\n\n" >> "step0_GxE_p_value_getters/${j}_${i}.sh"

        printf "python step0_TRACE_GxE_p_value_functions.py --i ${i} --j ${j}" >> "step0_GxE_p_value_getters/${j}_${i}.sh"
    done
done

for j in $(seq 0 4); do
    for i in $(seq 0 1999); do 
        bsub < step0_GxE_p_value_getters/${j}_${i}.sh
    done
done