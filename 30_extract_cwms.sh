#!/bin/bash
#$ -cwd
#$ -N extract_cwms
#$ -j y
#$ -pe smp 20
#$ -q main.q@badger.stanford.edu

# Source conda environment

source /storage/kaelanb/anaconda3/etc/profile.d/conda.sh

# activate chrombpnet since this one has the python installs

conda activate chrombpnet

set -u

FINEMO_DIR=/storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/hitcall/all_peaks
OUT_DIR=/storage/kaelanb/analysis/lucia/model/snail_modisco_test/all_peaks/hitcall/all_peaks

python extract_cwms.py --indir "${FINEMO_DIR}" --outdir "${OUT_DIR}"
