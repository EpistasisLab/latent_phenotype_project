#!/bin/bash

#commands=("python step0_compute_SV_p_values.py --chr 1 --pheno 13 --rsID rs11799887 --model PCA --count 3"
#          "python step0_compute_SV_p_values.py --chr 1 --pheno 13 --rsID rs11799887 --model PCA --count 4"
#          "python step0_compute_SV_p_values.py --chr 1 --pheno 13 --rsID rs11799887 --model PCA --count 5"
#          "python step0_compute_SV_p_values.py --chr 1 --pheno 13 --rsID rs11799887 --model PCA --count 6"
#          "python step0_compute_SV_p_values.py --chr 1 --pheno 13 --rsID rs11799887 --model PCA --count 7"
#          "python step0_compute_SV_p_values.py --chr 1 --pheno 13 --rsID rs11799887 --model PCA --count 8"
#          "python step0_compute_SV_p_values.py --chr 1 --pheno 13 --rsID rs11799887 --model PCA --count 9"
#          "python step0_compute_SV_p_values.py --chr 1 --pheno 13 --rsID rs11799887 --model PCA --count 10"
#          "python step0_compute_SV_p_values.py --chr 1 --pheno 13 --rsID rs11799887 --model PCA --count 11"
#          "python step0_compute_SV_p_values.py --chr 1 --pheno 13 --rsID rs11799887 --model PCA --count 12"
#)

#commands=("python step0_compute_SV_p_values.py --chr 3 --pheno 12 --rsID rs67183339 --model PCA --count 3"
#          "python step0_compute_SV_p_values.py --chr 3 --pheno 12 --rsID rs67183339 --model PCA --count 4"
#          "python step0_compute_SV_p_values.py --chr 3 --pheno 12 --rsID rs67183339 --model PCA --count 5"
#          "python step0_compute_SV_p_values.py --chr 3 --pheno 12 --rsID rs67183339 --model PCA --count 6"
#          "python step0_compute_SV_p_values.py --chr 3 --pheno 12 --rsID rs67183339 --model PCA --count 7"
#          "python step0_compute_SV_p_values.py --chr 3 --pheno 12 --rsID rs67183339 --model PCA --count 8"
#          "python step0_compute_SV_p_values.py --chr 3 --pheno 12 --rsID rs67183339 --model PCA --count 9"
#          "python step0_compute_SV_p_values.py --chr 3 --pheno 12 --rsID rs67183339 --model PCA --count 10"
#          "python step0_compute_SV_p_values.py --chr 3 --pheno 12 --rsID rs67183339 --model PCA --count 11"
#          "python step0_compute_SV_p_values.py --chr 3 --pheno 12 --rsID rs67183339 --model PCA --count 12"
#)

commands=("python step0_compute_SV_p_values.py --chr 14 --pheno 12 --rsID rs17486830 --model PCA --count 3"
          "python step0_compute_SV_p_values.py --chr 14 --pheno 12 --rsID rs17486830 --model PCA --count 4"
          "python step0_compute_SV_p_values.py --chr 14 --pheno 12 --rsID rs17486830 --model PCA --count 5"
          "python step0_compute_SV_p_values.py --chr 14 --pheno 12 --rsID rs17486830 --model PCA --count 6"
          "python step0_compute_SV_p_values.py --chr 14 --pheno 12 --rsID rs17486830 --model PCA --count 7"
          "python step0_compute_SV_p_values.py --chr 14 --pheno 12 --rsID rs17486830 --model PCA --count 8"
          "python step0_compute_SV_p_values.py --chr 14 --pheno 12 --rsID rs17486830 --model PCA --count 9"
          "python step0_compute_SV_p_values.py --chr 14 --pheno 12 --rsID rs17486830 --model PCA --count 10"
          "python step0_compute_SV_p_values.py --chr 14 --pheno 12 --rsID rs17486830 --model PCA --count 11"
          "python step0_compute_SV_p_values.py --chr 14 --pheno 12 --rsID rs17486830 --model PCA --count 12"
)

for cmd in "${commands[@]}"; do
    job_script="job_$(date +%Y%m%d%H%M%S%N).sh"
    cp step7g_sub_phenotype_analysis_template.sh $job_script  # Use cp command to copy the template
    sed -i "s|PLACEHOLDER|$cmd|g" $job_script
    if bsub < $job_script; then
        rm  $job_script
    else
        echo "Job submission failed for $job_script"
    fi
done