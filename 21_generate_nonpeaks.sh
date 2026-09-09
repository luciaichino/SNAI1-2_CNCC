#!/bin/bash
#$ -cwd
#$ -N chrombpnet_nonpeaks
#$ -pe smp 10

# Source conda environment

source /storage/kaelanb/anaconda3/etc/profile.d/conda.sh

# activate chrombpnet

conda activate chrombpnet

# Generate non-peak background regions

chrombpnet prep nonpeaks \
-g /storage/kaelanb/scripts/chrombpnet/hg38.fa \
-p /storage/kaelanb/analysis/lucia/data/data_forKaelan/20260607_cncc_training_atac_peaks.narrowpeak \
-c  /storage/kaelanb/scripts/chrombpnet/hg38.chrom.sizes \
-fl /storage/kaelanb/scripts/chrombpnet/fold_0.json \
-br /storage/kaelanb/scripts/chrombpnet/blacklist.bed.gz \
-o /storage/kaelanb/analysis/lucia/data/data_forKaelan/NONPEAK_20260607_cncc_training_atac_peaks
