#!/bin/bash
#BSUB -J step9d_get_QQ_plots
#BSUB -o step9d_get_QQ_plots.out
#BSUB -e step9d_get_QQ_plots.err

DIR="QQ_plots_getters_smoking"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for i in {0..15}
do
    printf '#!/bin/bash'"\n" > "QQ_plots_getters_smoking/get_QQ_plot_phenotype"$i".sh"
    printf '#BSUB -J step9d_get_QQ_plot_phenotype'$i"\n" >> "QQ_plots_getters_smoking/get_QQ_plot_phenotype"$i".sh"
    printf '#BSUB -o QQ_plots_getters_smoking/get_QQ_plot_phenotype'$i".out\n" >> "QQ_plots_getters_smoking/get_QQ_plot_phenotype"$i".sh"
    printf '#BSUB -e QQ_plots_getters_smoking/get_QQ_plot_phenotype'$i".err\n" >> "QQ_plots_getters_smoking/get_QQ_plot_phenotype"$i".sh"
    printf '#BSUB -R "rusage[mem=30000MB]"'"\n" >> "QQ_plots_getters_smoking/get_QQ_plot_phenotype"$i".sh"
    printf '#BSUB -M 30000MB'"\n" >> "QQ_plots_getters_smoking/get_QQ_plot_phenotype"$i".sh"
    printf 'source activate torch_env2'"\n\n" >> "QQ_plots_getters_smoking/get_QQ_plot_phenotype"$i".sh"
    printf 'python step0_QQ_plot_getter.py --env smoking --pheno '$i >> "QQ_plots_getters_smoking/get_QQ_plot_phenotype"$i".sh"
done

for i in {0..15}
do 
    bsub < QQ_plots_getters_smoking/get_QQ_plot_phenotype"$i".sh
done