#!/bin/bash
#$ -cwd
#$ -N convert_h5_npz
#$ -pe smp 10
#$ -q main.q@badger.stanford.edu

# Source conda environment

source /storage/kaelanb/anaconda3/etc/profile.d/conda.sh

# activate h5_to_npz: this is a temp conda environment for the conversion only

conda activate h5_to_npz

# Run the conversion

python convert_chrombpnet_h5_to_npz.py \
    --input /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/bw/all_peaks.counts_scores.h5 \
    --ohe_output /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/bw/all_peaks.counts_scores_ohe.npz \
    --attr_output /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/bw/all_peaks.counts_scores_attr.npz

python convert_chrombpnet_h5_to_npz.py \
    --input /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/bw/all_peaks.profile_scores.h5 \
    --ohe_output /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/bw/all_peaks.profile_scores_ohe.npz \
    --attr_output /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/bw/all_peaks.profile_scores_attr.npz
