#!/bin/bash
#$ -cwd
#$ -N modisco
#$ -pe smp 60
#$ -j y
#$ -q main.q@raccoon.stanford.edu
#$ -o relaxed.log

# Source conda environment

source /storage/kaelanb/anaconda3/etc/profile.d/conda.sh

# activate modisco lite environment

conda activate tfmodisco-lite

python run_modisco.py \
    -s /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/bw/all_peaks.counts_scores_ohe.npz \
    -a /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/bw/all_peaks.counts_scores_attr.npz \
    -n 500000 \
    -o /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/modisco/counts/relaxed/relaxed_fdr005.h5 \
    --fdr 0.05

modisco report \
    -i /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/modisco/counts/relaxed/relaxed_fdr005.h5 \
    -o /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/modisco/counts/relaxed/ \
    -s /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/modisco/counts/relaxed/
sed -i 's|/storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/modisco/counts/relaxed/|./|g' \
    /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/modisco/counts/relaxed/motifs.html
