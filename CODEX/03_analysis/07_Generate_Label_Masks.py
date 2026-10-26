import argparse
import os

import numpy as np
import pandas as pd
import tifffile

parser = argparse.ArgumentParser()
parser.add_argument("--segmentation", required=True, help="Mesmer segmentation_cell.tif")
parser.add_argument("--labels", required=True, help="<sample>_labels.csv")
parser.add_argument("--column", default="cell_type_code",
                    choices=["cell_type_code", "cn_code"])
parser.add_argument("--output", required=True, help="output prefix (.npy and .tif)")
args = parser.parse_args()

# ==============================================================================
# 1. Load segmentation and labels
# ==============================================================================
mask = tifffile.imread(args.segmentation).astype(np.int64)
labels = pd.read_csv(args.labels).dropna(subset=[args.column])
labels = labels[labels["CellID"] <= mask.max()]

# ==============================================================================
# 2. Map cell IDs to label codes with a lookup table
# ==============================================================================
lut = np.zeros(mask.max() + 1, dtype=np.uint8)
lut[labels["CellID"].to_numpy(dtype=np.int64)] = labels[args.column].to_numpy(dtype=np.uint8)
out = lut[mask]

# ==============================================================================
# 3. Save
# ==============================================================================
os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
np.save(args.output + ".npy", out)
tifffile.imwrite(args.output + ".tif", out)
print(f"Saved {args.output}.npy/.tif: {labels.shape[0]} labelled cells")
