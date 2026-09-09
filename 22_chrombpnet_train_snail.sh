#!/bin/bash

# Source conda environment

source ~/anaconda3/etc/profile.d/conda.sh

# activate chrombpnet

conda activate chrombpnet

# Specify which GPU to use (e.g., GPU 0, or GPU 1, etc.)
export CUDA_VISIBLE_DEVICES=1  # Use GPU 0

# Train single task ChromBPNet models

chrombpnet pipeline \
    -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_DMSO_merged.bam \
    -d "ATAC" \
    -g ~/analysis/genomes/hg38.fa \
    -c ~/analysis/genomes/hg38.chrom.sizes \
    -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
    -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
    -fl ~/analysis/genomes/fold_0.json \
    -b ~/analysis/genomes/ENCSR868FGK_bias_fold_0.h5 \
    -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks

chrombpnet pipeline \
    -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_3h_merged.bam \
    -d "ATAC" \
    -g ~/analysis/genomes/hg38.fa \
    -c ~/analysis/genomes/hg38.chrom.sizes \
    -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
    -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
    -fl ~/analysis/genomes/fold_0.json \
    -b ~/analysis/genomes/ENCSR868FGK_bias_fold_0.h5 \
    -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/dtag_3h/all_peaks

# Generate contrib

chrombpnet contribs_bw \
    -m ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/models/chrombpnet_nobias.h5 \
    -r ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
    -g ~/analysis/genomes/hg38.fa \
    -c ~/analysis/genomes/hg38.chrom.sizes \
    -op ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/bw/all_peaks

chrombpnet contribs_bw \
    -m ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/dtag_3h/all_peaks/models/chrombpnet_nobias.h5 \
    -r ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
    -g ~/analysis/genomes/hg38.fa \
    -c ~/analysis/genomes/hg38.chrom.sizes \
    -pc counts \
    -op ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/dtag_3h/all_peaks/bw/all_peaks

# Pred bw

chrombpnet pred_bw \
    -bm ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/models/bias_model_scaled.h5 \
    -cm ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/models/chrombpnet.h5 \
    -cmb ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/models/chrombpnet_nobias.h5 \
    -r ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
    -g ~/analysis/genomes/hg38.fa \
    -c ~/analysis/genomes/hg38.chrom.sizes \
    -op ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/bw/all_peaks

chrombpnet pred_bw \
    -bm ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/dtag_3h/all_peaks/models/bias_model_scaled.h5 \
    -cm ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/dtag_3h/all_peaks/models/chrombpnet.h5 \
    -cmb ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/dtag_3h/all_peaks/models/chrombpnet_nobias.h5 \
    -r ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
    -g ~/analysis/genomes/hg38.fa \
    -c ~/analysis/genomes/hg38.chrom.sizes \
    -op ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/all_peaks/dtag_3h/all_peaks/bw/all_peaks
