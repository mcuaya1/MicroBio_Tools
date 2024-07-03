#!/bin/bash

CONDA_PATH=~/miniconda3/etc/profile.d/conda.sh
MAIN_QZA=/home/mario/Desktop/working-dir/Fungal/FUN_L6.qza
OUTPUTDIR=/home/mario/Desktop/working-dir/Fungal/ROOTS/ALL
GROUPED="ALL_GROUPED"
CATEGORY="roots"
MAP_FILE=/home/mario/Desktop/working-dir/Fungal/fungal-master-map.tsv
FILTER_STRING="'T1Tm0_CAR', 'T2Tm0_CAR', 'T3Tm0_CAR', 'T4Tm0_CAR', 'T5Tm0_CAR', 'T6Tm0_CAR', 'T7Tm0_CAR', 'T8Tm0_CAR', 'T1Tm0_C35', 'T2Tm0_C35', 'T3Tm0_C35', 'T4Tm0_C35', 'T5Tm0_C35', 'T6Tm0_C35', 'T7Tm0_C35', 'T8Tm0_C35', 'T1Tm154_CAR', 'T2Tm154_CAR', 'T3Tm154_CAR', 'T4Tm154_CAR', 'T5Tm154_CAR', 'T6Tm154_CAR', 'T7Tm154_CAR', 'T8Tm154_CAR', 'T1Tm154_C35', 'T2Tm154_C35', 'T3Tm154_C35', 'T4Tm154_C35', 'T5Tm154_C35', 'T6Tm154_C35', 'T7Tm154_C35', 'T8Tm154_C35', 'T9Tm0', 'T10Tm0'"

source ${CONDA_PATH}
#Used to group samples with the same metadata tag in a specific colum()
qiime feature-table filter-samples \
--i-table FUN_L6.qza \
--m-metadata-file fungal-master-map.tsv \
--p-where "PHASE_3_Liquid_Greenhouse_SOIL IN (${FILTER_STRING})" \
--o-filtered-table ${OUTPUTDIR}/treatment_filtered_${CATEGORY}_${GROUPED}.qza

