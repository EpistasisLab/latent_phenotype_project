#!/bin/bash
#BSUB -J step3d_run_king_subjobs
#BSUB -o step3d_run_king_subjobs.out
#BSUB -e step3d_run_king_subjobs.error

bsub < king_subjobs/get_kinships_1.sh
bsub < king_subjobs/get_kinships_2.sh
bsub < king_subjobs/get_kinships_3.sh
bsub < king_subjobs/get_kinships_4.sh
bsub < king_subjobs/get_kinships_5.sh
bsub < king_subjobs/get_kinships_6.sh
bsub < king_subjobs/get_kinships_7.sh
bsub < king_subjobs/get_kinships_8.sh
bsub < king_subjobs/get_kinships_9.sh
bsub < king_subjobs/get_kinships_10.sh

bsub < king_subjobs/get_kinships_12.sh
bsub < king_subjobs/get_kinships_13.sh
bsub < king_subjobs/get_kinships_14.sh
bsub < king_subjobs/get_kinships_15.sh
bsub < king_subjobs/get_kinships_16.sh
bsub < king_subjobs/get_kinships_17.sh
bsub < king_subjobs/get_kinships_18.sh
bsub < king_subjobs/get_kinships_19.sh
bsub < king_subjobs/get_kinships_110.sh

bsub < king_subjobs/get_kinships_23.sh
bsub < king_subjobs/get_kinships_24.sh
bsub < king_subjobs/get_kinships_25.sh
bsub < king_subjobs/get_kinships_26.sh
bsub < king_subjobs/get_kinships_27.sh
bsub < king_subjobs/get_kinships_28.sh
bsub < king_subjobs/get_kinships_29.sh
bsub < king_subjobs/get_kinships_210.sh

bsub < king_subjobs/get_kinships_34.sh
bsub < king_subjobs/get_kinships_35.sh
bsub < king_subjobs/get_kinships_36.sh
bsub < king_subjobs/get_kinships_37.sh
bsub < king_subjobs/get_kinships_38.sh
bsub < king_subjobs/get_kinships_39.sh
bsub < king_subjobs/get_kinships_310.sh

bsub < king_subjobs/get_kinships_45.sh
bsub < king_subjobs/get_kinships_46.sh
bsub < king_subjobs/get_kinships_47.sh
bsub < king_subjobs/get_kinships_48.sh
bsub < king_subjobs/get_kinships_49.sh
bsub < king_subjobs/get_kinships_410.sh

bsub < king_subjobs/get_kinships_56.sh
bsub < king_subjobs/get_kinships_57.sh
bsub < king_subjobs/get_kinships_58.sh
bsub < king_subjobs/get_kinships_59.sh
bsub < king_subjobs/get_kinships_510.sh

bsub < king_subjobs/get_kinships_67.sh
bsub < king_subjobs/get_kinships_68.sh
bsub < king_subjobs/get_kinships_69.sh
bsub < king_subjobs/get_kinships_610.sh

bsub < king_subjobs/get_kinships_78.sh
bsub < king_subjobs/get_kinships_79.sh
bsub < king_subjobs/get_kinships_710.sh

bsub < king_subjobs/get_kinships_89.sh
bsub < king_subjobs/get_kinships_810.sh

bsub < king_subjobs/get_kinships_910.sh

