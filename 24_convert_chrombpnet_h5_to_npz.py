import hdf5plugin
import h5py
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Convert .h5 to .npz for TF-MoDISco-lite input.")
parser.add_argument('--input', required=True, help='Path to .h5 file')
parser.add_argument('--ohe_output', required=True, help='Path to output one-hot encoding .npz')
parser.add_argument('--attr_output', required=True, help='Path to output attribution scores .npz')

args = parser.parse_args()

# Open the .h5 file
with h5py.File(args.input, "r") as f:
    # Extract the one-hot encoding (raw) for sequences
    one_hot = f["/raw/seq"][:]

    # Extract the hypothetical contribution scores (shap) for attribution
    attribution_scores = f["/shap/seq"][:]  # This is the DeepSHAP contribution scores

# Save them as separate .npz files
np.savez(args.ohe_output, one_hot)
np.savez(args.attr_output, attribution_scores)

print(f"✅ Successfully saved one-hot encoding to {args.ohe_output}")
print(f"✅ Successfully saved attribution scores to {args.attr_output}")
