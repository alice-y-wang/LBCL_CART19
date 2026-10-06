import argparse
import os

import tifffile
from wsireg.wsireg2d import WsiReg2D

CODEX_RESOLUTION = 0.5069   # um / pixel
HE_RESOLUTION = 0.2485      # um / pixel

parser = argparse.ArgumentParser()
parser.add_argument("--name", required=True, help="sample / slide name")
parser.add_argument("--codex", required=True, help="CODEX .qptiff or .tif")
parser.add_argument("--he", required=True, help="H&E .qptiff or .tif")
parser.add_argument("--outdir", required=True)
args = parser.parse_args()
os.makedirs(args.outdir, exist_ok=True)

# ==============================================================================
# 1. Convert .qptiff to .tif
# ==============================================================================
def to_tif(path, out_path):
    if not path.endswith(".qptiff"):
        return path
    with tifffile.TiffFile(path) as tif:
        img = tif.series[0].asarray()
    tifffile.imwrite(out_path, img, bigtiff=True)
    return out_path

codex_tif = to_tif(args.codex, os.path.join(args.outdir, f"{args.name}_CODEX.tif"))
he_tif = to_tif(args.he, os.path.join(args.outdir, f"{args.name}_HE.tif"))

# ==============================================================================
# 2. Register H&E to CODEX
# ==============================================================================
reg_graph = WsiReg2D(args.name, os.path.join(args.outdir, "Output"))

reg_graph.add_modality(
    "CODEX", codex_tif, image_res=CODEX_RESOLUTION,
    preprocessing={"image_type": "FL", "ch_indices": [0],
                   "as_uint8": True, "contrast_enhance": True},
)
reg_graph.add_modality(
    "HE", he_tif, image_res=HE_RESOLUTION,
    preprocessing={"image_type": "BF", "as_uint8": True, "invert_intensity": True},
)
reg_graph.add_reg_path("HE", "CODEX", thru_modality=None, reg_params=["rigid", "affine"])

reg_graph.register_images()
reg_graph.save_transformations()
reg_graph.transform_images(file_writer="ome.tiff")
print("Registration done")
