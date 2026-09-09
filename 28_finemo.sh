#!/bin/bash
#$ -cwd
#$ -N finemo
#$ -pe smp 20
#$ -q main.q@badger.stanford.edu

# Source conda environment

source /storage/kaelanb/anaconda3/etc/profile.d/conda.sh

# activate chrombpnet

conda activate finemo

# Extract sequences and contributions from tfmodisco-lite input .npz files: the same used for TFmodisco-lite.

finemo extract-regions-modisco-fmt \
-s /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/bw/all_peaks.counts_scores_ohe.npz \
-a /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/bw/all_peaks.counts_scores_attr.npz \
-o /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/hitcall/all_peaks \
-p /storage/kaelanb/analysis/lucia/data/data_forKaelan/20260607_cncc_training_atac_peaks.narrowpeak

# Call hits

finemo call-hits \
-r /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/hitcall/all_peaks.npz \
-m /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/modisco/counts/relaxed/relaxed_fdr01.h5 \
-o /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/hitcall/all_peaks

# Generate report

finemo report \
-r /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/hitcall/all_peaks.npz \
-H /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/hitcall/all_peaks/ \
-o /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/hitcall/all_peaks/ \
-m /storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/modisco/counts/relaxed/relaxed_fdr01.h5
