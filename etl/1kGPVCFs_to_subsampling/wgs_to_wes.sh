#!/bin/bash --login

# How to run
# bash -i script.sh /!\

#############
# Ressources
# https://github.com/conda/conda/issues/13002
# https://github.com/conda/conda/issues/7126

#############
# Error handling
set -e

#############
# Variables
chromosome=$1
sourcedir="/data/1KG_GRCh38/"
wgsvcf=$sourcedir"/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chromosome}.recalibrated_variants.vcf.gz"
bed=$2

#############
# Activating conda environement
conda activate ${CONDA_PREFIX}/envs/labaz

#############
# Extracting exome
echo "$0:L${LINENO}:Extracting exome from chromosome $chromosome vcf."
{
if [ ! -f $wgsvcf ]; then
    echo "$0:L${LINENO}:File not found."; echo $wgsvcf
    exit 0
fi
}
wesvcf=$(echo $wgsvcf | sed -e 's,20201028_CCDG,exons_geneitx_20201028_CCDG,' -e 's,.vcf.gz,.vcf,')
echo "$0:L${LINENO}:Writing to:" $wesvcf
bcftools view -R ${bed} $wgsvcf | bcftools sort -Oz - > ${wesvcf}.gz
echo "$0:L${LINENO}:Script finished."