#!/bin/bash
#BSUB -J step9c_get_significant_QTLs
#BSUB -o step9c_get_significant_QTLs.out
#BSUB -e step9c_get_significant_QTLs.err

DIR="binary_HF_QTL_getters"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for i in {1..22} "MT" "X" "XY" "Y"
do
   for j in "cong_HF" "CMyo" "LV_HF" "Unknown_HF" "any_HF"
   do
      for k in "0" "2000" "4000" "6000" "8000" "10000" "12000" "14000"
      do 
         printf '#!/bin/bash'"\n" > "binary_HF_QTL_getters/get_chr"$i"_"$j"_jstart"$k"_binary_HF_QTLs.sh"
         printf '#BSUB -J get_chr'$i'_pheno'$j"_jstart"$k"_binary_HF_QTLs\n" >> "binary_HF_QTL_getters/get_chr"$i"_"$j"_jstart"$k"_binary_HF_QTLs.sh"
         printf '#BSUB -o binary_HF_QTL_getters/get_chr'$i'_pheno'$j"_jstart"$k"_binary_HF_QTLs.out\n" >> "binary_HF_QTL_getters/get_chr"$i"_"$j"_jstart"$k"_binary_HF_QTLs.sh"
         printf '#BSUB -e binary_HF_QTL_getters/get_chr'$i'_pheno'$j"_jstart"$k"_binary_HF_QTLs.err\n" >> "binary_HF_QTL_getters/get_chr"$i"_"$j"_jstart"$k"_binary_HF_QTLs.sh"
         printf '#BSUB -R "rusage[mem=15000MB]"\n' >> "binary_HF_QTL_getters/get_chr"$i"_"$j"_jstart"$k"_binary_HF_QTLs.sh"
         printf '#BSUB -M 15000MB\n' >> "binary_HF_QTL_getters/get_chr"$i"_"$j"_jstart"$k"_binary_HF_QTLs.sh"
         printf 'source activate torch_env2'"\n\n" >> "binary_HF_QTL_getters/get_chr"$i"_"$j"_jstart"$k"_binary_HF_QTLs.sh"
         printf 'python step0_binary_HF_QTL_getter.py --chr '$i" --pheno "$j" --jstart "$k"\n">> "binary_HF_QTL_getters/get_chr"$i"_"$j"_jstart"$k"_binary_HF_QTLs.sh"
      done
   done
done

for i in {1..22} "MT" "X" "XY" "Y"
do
   for j in "cong_HF" "CMyo" "LV_HF" "Unknown_HF" "any_HF"
   do
      for k in "0" "2000" "4000" "6000" "8000" "10000" "12000" "14000"
      do
         bsub < binary_HF_QTL_getters/get_chr"$i"_"$j"_jstart"$k"_binary_HF_QTLs.sh
      done
   done
done
