#!/usr/bin/env python
"""
run_modisco_relaxed.py

A drop-in replacement for `modisco motifs` that exposes target_seqlet_fdr,
which is hardcoded (0.05) in the tfmodisco-lite CLI and not otherwise
changeable without editing the package source.

Does NOT modify your original .npz files, tfmodisco-lite package, or
environment in any way -- it just imports the already-installed modiscolite
package and calls its Python API with a different fdr value.

USAGE (mirrors your existing `modisco motifs` command):

    python run_modisco_relaxed.py \
        -s all_peaks_no_prom_w_10perc.counts_scores_ohe.npz \
        -a all_peaks_no_prom_w_10perc.counts_scores_attr.npz \
        -n 500000 \
        -o modisco_results_fdr02.h5 \
        --fdr 0.2

Then run `modisco report` on the resulting .h5 exactly as you did before.

After this finishes, look at the printed seqlet counts:
  - If the negative metacluster is still empty/tiny even at a loose FDR
    (e.g. --fdr 0.3), that's evidence this is a real-signal issue, not a
    threshold issue.
  - If a negative metacluster appears once you loosen the FDR, that
    confirms it was a threshold problem -- use whichever fdr value gave
    you a usable, stable set of patterns.
"""

import argparse
import sys
import numpy as np


def load_array(npz_path, label):
    data = np.load(npz_path)
    key = data.files[0]
    arr = data[key]
    print(f"[{label}] loaded '{npz_path}' -> array '{key}' with shape {arr.shape}, dtype {arr.dtype}")
    return arr


def main():
    parser = argparse.ArgumentParser(
        description="Run tfmodisco-lite with a custom target_seqlet_fdr (not exposed in the normal CLI)."
    )
    parser.add_argument("-s", "--sequences", required=True, help="Path to one-hot encoded sequences .npz (same file you pass to `modisco motifs -s`)")
    parser.add_argument("-a", "--attributions", required=True, help="Path to hypothetical contribution scores .npz (same file you pass to `modisco motifs -a`)")
    parser.add_argument("-n", "--max_seqlets", type=int, default=500000, help="Max seqlets per metacluster (same as `modisco motifs -n`). Default: 500000")
    parser.add_argument("-o", "--output", required=True, help="Output .h5 path (same format as `modisco motifs -o`, usable directly with `modisco report`)")
    parser.add_argument("-w", "--window", type=int, default=20, help="Sliding window size used for seqlet calling. Default: 20")
    parser.add_argument("--flank_size", type=int, default=5, help="Flank size added around each seqlet. Default: 5")
    parser.add_argument("--n_leiden", type=int, default=2, help="Number of Leiden clustering runs (same as `modisco motifs -l`). Default: 2")
    parser.add_argument("--fdr", type=float, default=0.2, help="target_seqlet_fdr. Default in tfmodisco-lite CLI is hardcoded to 0.05. Raise this (e.g. 0.15-0.3) to allow MORE seqlets through -- try this if you're missing a negative metacluster. Lowering it gives FEWER, stricter seqlets.")
    args = parser.parse_args()

    try:
        import modiscolite
        from modiscolite.tfmodisco import TFMoDISco
        from modiscolite.io import save_hdf5
    except ImportError:
        print("ERROR: could not import modiscolite. Make sure you're running this in the same "
              "environment/conda env where `modisco motifs` normally works (i.e. wherever "
              "tfmodisco-lite is installed). This script does not install or change anything.",
              file=sys.stderr)
        sys.exit(1)

    print(f"\nUsing tfmodisco-lite version: {getattr(modiscolite, '__version__', 'unknown')}\n")

    one_hot = load_array(args.sequences, "one_hot")
    hyp_scores = load_array(args.attributions, "hypothetical_contribs")

    print("Original shapes:")
    print("one_hot:", one_hot.shape)
    print("hyp_scores:", hyp_scores.shape)

# Convert from (N, 4, L) -> (N, L, 4), as expected by TFMoDISco API
    one_hot = np.transpose(one_hot, (0, 2, 1)).astype(np.float32)
    hyp_scores = np.transpose(hyp_scores, (0, 2, 1)).astype(np.float32)

    print("After transpose:")
    print("one_hot:", one_hot.shape)
    print("hyp_scores:", hyp_scores.shape)

    # Sanity check: shapes should match on the sequence-count and length axes
    if one_hot.shape != hyp_scores.shape:
        print(f"WARNING: one_hot shape {one_hot.shape} does not match hypothetical_contribs shape "
              f"{hyp_scores.shape}. These are normally expected to match exactly "
              f"(same N, same 4-letter axis, same length axis). Double check you passed the "
              f"correct pair of files before proceeding.", file=sys.stderr)

    print(f"\nRunning TFMoDISco with target_seqlet_fdr={args.fdr}, "
          f"max_seqlets_per_metacluster={args.max_seqlets}, sliding_window_size={args.window}, "
          f"flank_size={args.flank_size}, n_leiden_runs={args.n_leiden}\n")

    pos_patterns, neg_patterns = TFMoDISco(
        hypothetical_contribs=hyp_scores,
        one_hot=one_hot,
        max_seqlets_per_metacluster=args.max_seqlets,
        sliding_window_size=args.window,
        flank_size=args.flank_size,
        target_seqlet_fdr=args.fdr,
        n_leiden_runs=args.n_leiden,
        verbose=True,
    )

    n_pos_patterns = len(pos_patterns) if pos_patterns is not None else 0
    n_neg_patterns = len(neg_patterns) if neg_patterns is not None else 0

    def total_seqlets(patterns):
        if not patterns:
            return 0
        total = 0
        for p in patterns:
            try:
                total += len(p.seqlets)
            except AttributeError:
                pass
        return total

    print("\n===== RESULTS SUMMARY =====")
    print(f"Positive patterns found: {n_pos_patterns} (total seqlets: {total_seqlets(pos_patterns)})")
    print(f"Negative patterns found: {n_neg_patterns} (total seqlets: {total_seqlets(neg_patterns)})")
    print("============================\n")

    if n_neg_patterns == 0:
        print("Still no negative patterns at this FDR. Try rerunning with a higher --fdr "
              "(e.g. 0.3), and/or check the raw attribution score distribution (min/percentiles) "
              "to see if there's meaningful negative signal in this model's attributions at all.")

    save_hdf5(args.output, pos_patterns, neg_patterns, args.window)
    print(f"Saved results to: {args.output}")
    print(f"Next step: modisco report -i {args.output} -o <report_dir> -s <report_dir>")


if __name__ == "__main__":
    main()
