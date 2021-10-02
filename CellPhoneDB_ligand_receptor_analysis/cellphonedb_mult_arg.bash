#!/bin/bash

SAMPLE=$1
granularity=$2

echo "**************************"
echo "Running CellPhoneDB  $(date)"
echo "Patient:  ${SAMPLE}"
echo "**************************"
cd data/CellPhoneDB/${SAMPLE}
cellphonedb method statistical_analysis ${SAMPLE}_${granularity}_meta.txt ${SAMPLE}_count.txt --output-path cell_type_${granularity} --counts-data gene_name --threads 4
#heatmap
cellphonedb plot heatmap_plot ${SAMPLE}_${granularity}_meta.txt --pvalues-path cell_type_${granularity}/pvalues.txt --output-path cell_type_${granularity} --count-name heatmap_${SAMPLE}_${granularity}.pdf
#plotting
cellphonedb plot dot_plot --means-path cell_type_${granularity}/means.txt --pvalues-path cell_type_${granularity}/pvalues.txt --output-path cell_type_${granularity} --output-name dotplot_${SAMPLE}_${granularity}.pdf
echo "**************************"
echo "Finished CellPhoneDB  $(date)"
