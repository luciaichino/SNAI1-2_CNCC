#!/usr/bin/env bash
# Run the whole pipeline. Each step can also be run on its own after sourcing config.sh.
set -euo pipefail
cd "$(dirname "$0")"

# The config.sh file contains the paths to the input files and output directories. 
# You need to edit it before running this script.
source ./config.sh

bash    scripts/1_define_elements.sh
bash    scripts/2_covariates.sh
Rscript scripts/3_match_atac_background.R
Rscript scripts/4_match_kmer_background.R
bash    scripts/5_liftover_hg19.sh
bash    scripts/6_score.sh
Rscript scripts/7_statistics.R
Rscript scripts/8_figure.R

# The output PDF should be in $FIG/fig_final_comparison.pdf
