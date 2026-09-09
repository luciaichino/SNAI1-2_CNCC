#!/usr/bin/env python
"""
Kaelan Brennan
Extract Fi-NeMo CWMs (motif_cwms.npy + motif_data.tsv) into a tidy TSV.

IMPORTANT: motif_cwms.npy stores each motif's CWM on an internal scale used
for Fi-NeMo's shape-based (cosine-similarity) hit calling -- NOT necessarily
a scale that's comparable across motifs in absolute amplitude. motif_data.tsv's
`motif_scale` column is "the motif scaling factor, used to normalize by motif
importance" (per the Fi-NeMo README), so to recover real, cross-motif-comparable
contribution amplitudes, multiply each motif's raw CWM by its own motif_scale.
This script outputs both the raw and motif_scale-corrected values so you can
compare them directly against independently-computed importance metrics
(e.g. summed real contribution scores at each hit's true coordinates).

Usage:
    python extract_cwms.py --indir /path/to/finemo/outdir --outdir /path/to/out
"""

import argparse
import os

import numpy as np
import pandas as pd

BASES = ["A", "C", "G", "T"]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--indir", required=True, help="Fi-NeMo output directory "
                   "(contains motif_cwms.npy and motif_data.tsv)")
    p.add_argument("--outdir", required=True, help="Directory to write TSVs to")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    cwms = np.load(os.path.join(args.indir, "motif_cwms.npy"))  # (n, 4, w)
    motif_data = pd.read_csv(os.path.join(args.indir, "motif_data.tsv"), sep="\t")

    # Array row order corresponds to ascending motif_id
    motif_data = motif_data.sort_values("motif_id").reset_index(drop=True)

    if cwms.shape[0] != len(motif_data):
        raise ValueError(
            f"Row count mismatch: motif_cwms.npy has {cwms.shape[0]} motifs, "
            f"motif_data.tsv has {len(motif_data)} rows."
        )

    full_records = []
    trimmed_records = []

    for i, row in motif_data.iterrows():
        mat = cwms[i]  # (4, w)
        width = mat.shape[1]

        motif_id = row["motif_id"]
        motif_name = row["motif_name"]
        strand = row.get("strand", None)
        m_start = int(row["motif_start"]) if "motif_start" in row else 0
        m_end = int(row["motif_end"]) if "motif_end" in row else width
        motif_scale = float(row["motif_scale"]) if "motif_scale" in row else 1.0

        for pos in range(width):
            rec = {
                "motif_id": motif_id,
                "motif_name": motif_name,
                "strand": strand,
                "position": pos,
                "motif_scale": motif_scale,
            }
            for bi, b in enumerate(BASES):
                raw_val = mat[bi, pos]
                rec[b] = raw_val                      # raw npy value
                rec[f"{b}_scaled"] = raw_val * motif_scale  # motif_scale-corrected
            full_records.append(rec)

            if m_start <= pos < m_end:
                trec = dict(rec)
                trec["position"] = pos - m_start  # re-zero to trimmed window
                trimmed_records.append(trec)

    full_df = pd.DataFrame(full_records)
    trimmed_df = pd.DataFrame(trimmed_records)

    full_path = os.path.join(args.outdir, "cwms_long_full.tsv")
    trimmed_path = os.path.join(args.outdir, "cwms_long_trimmed.tsv")

    full_df.to_csv(full_path, sep="\t", index=False)
    trimmed_df.to_csv(trimmed_path, sep="\t", index=False)

    print(f"Wrote {len(full_df)} rows ({motif_data.shape[0]} motifs) to {full_path}")
    print(f"Wrote {len(trimmed_df)} rows to {trimmed_path}")


if __name__ == "__main__":
    main()
