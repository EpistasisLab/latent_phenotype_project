#!/bin/bash
#BSUB -J step9f_get_QQ_plots_normal_GWAS
#BSUB -o step9f_get_QQ_plots_normal_GWAS.out
#BSUB -e step9f_get_QQ_plots_normal_GWAS.err

DIR="QQ_plots_getters_normal_GWAS"
if test ! -d $DIR; then
    mkdir "$DIR"
fi

for i in any_HF LV_HF Unknown_HF CMyo cong_HF
do
    printf '#!/bin/bash'"\n" > "QQ_plots_getters_normal_GWAS/get_QQ_plot_phenotype"$i".sh"
    printf '#BSUB -J step9d_get_QQ_plot_phenotype'$i"\n" >> "QQ_plots_getters_normal_GWAS/get_QQ_plot_phenotype"$i".sh"
    printf '#BSUB -o QQ_plots_getters_normal_GWAS/get_QQ_plot_phenotype'$i".out\n" >> "QQ_plots_getters_smoking/get_QQ_plot_phenotype"$i".sh"
    printf '#BSUB -e QQ_plots_getters_normal_GWAS/get_QQ_plot_phenotype'$i".err\n" >> "QQ_plots_getters_smoking/get_QQ_plot_phenotype"$i".sh"
    printf '#BSUB -R "rusage[mem=30000MB]"'"\n" >> "QQ_plots_getters_normal_GWAS/get_QQ_plot_phenotype"$i".sh"
    printf '#BSUB -M 30000MB'"\n" >> "QQ_plots_getters_normal_GWAS/get_QQ_plot_phenotype"$i".sh"
    printf 'source activate torch_env2'"\n\n" >> "QQ_plots_getters_normal_GWAS/get_QQ_plot_phenotype"$i".sh"
    printf 'python step0_QQ_plot_getter.py --env normal_GWAS --pheno '$i' --GWAS normal' >> "QQ_plots_getters_normal_GWAS/get_QQ_plot_phenotype"$i".sh"
done

for i in any_HF LV_HF Unknown_HF CMyo cong_HF
do 
    bsub < QQ_plots_getters_normal_GWAS/get_QQ_plot_phenotype"$i".sh
done