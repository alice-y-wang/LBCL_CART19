import argparse

import dask.array as da
import napari
import numpy as np
import pandas as pd
import tifffile
import zarr

# ==============================================================================
# 1. Label colors
# ==============================================================================
CELLTYPE_COLORS = {
    1: "#D9D9D9",   # Lymphoma
    2: "#F9D448",   # Macrophage
    3: "#654321",   # CD14+TIMs
    4: "#FF9007",   # CD14+CD16+TIMs
    5: "#E94310",   # cDC
    6: "#E7298A",   # Granulocytes
    7: "#984EA3",   # T-regulatory
    8: "#1F9CDA",   # CD4 Naive-like
    9: "#112EB0",   # CD4 TEM
    10: "#917EE6",  # CD4T effector-exhausted
    11: "#5BAE03",  # CD8 TEM
    12: "#6FF588",  # CD8T effector-exhausted
    13: "#808080",  # Proliferating T
    14: "#00D6B6",  # Endothelial
    15: "#CFC691",  # Stromal
}
CN_COLORS = {
    1: "#FB8072", 2: "#8DD3C7", 3: "#BC80BD", 4: "#D9D9D9", 5: "#BEBADA",
    6: "#80B1D3", 7: "#FF7F00", 8: "#FDB462", 9: "#E7298A", 10: "#969696",
}


def highlight(colors, keep_codes, grey="#D9D9D9"):
    """Color only selected cell types (e.g. [11, 12] for CD8+ T cells)."""
    return {k: (v if k in keep_codes else grey) for k, v in colors.items()}

# ==============================================================================
# 2. Helpers
# ==============================================================================
def to_label_pyramid(mask, n_levels=4, downscale=2, chunks=2048):
    """Nearest-neighbour multiscale pyramid for a categorical mask."""
    mask = mask.astype(np.uint8 if mask.max() < 256 else np.uint16, copy=False)
    pyramid = [da.from_array(mask, chunks=chunks)]
    for _ in range(n_levels - 1):
        pyramid.append(pyramid[-1][::downscale, ::downscale])
    return pyramid


def lazy_view(tiff_path, channel_names=None, viewer=None, labs=None, cols=None,
              cn_labs=None, cn_cols=None, he_path=None, background_white=True):
    """Load a multiscale CODEX OME-TIFF lazily and add mask / H&E layers."""
    with tifffile.TiffFile(tiff_path, is_ome=False) as tiff:
        n_levels = len(tiff.series[0].levels)
        channel_axis = 0 if len(tiff.series[0].shape) == 3 else None

    z = zarr.open(tifffile.imread(tiff_path, aszarr=True), mode="r")
    pyramid = [da.from_zarr(z[str(i)]) for i in range(n_levels)]

    viewer = viewer if viewer is not None else napari.Viewer()
    viewer.add_image(pyramid, multiscale=True, channel_axis=channel_axis,
                     visible=False, name=channel_names, blending="additive")
    viewer.scale_bar.visible = True
    viewer.scale_bar.unit = "um"

    if labs is not None:
        viewer.add_labels(to_label_pyramid(labs), multiscale=True,
                          name="Cell type mask", colormap=cols)
    if cn_labs is not None:
        viewer.add_labels(to_label_pyramid(cn_labs), multiscale=True,
                          name="Neighborhood mask", colormap=cn_cols)
    if he_path is not None:
        viewer.add_image(tifffile.imread(he_path), name="H&E", blending="additive")
    if background_white:
        viewer.theme = "light"
    return viewer


# ==============================================================================
# 3. Launch viewer
# ==============================================================================
    parser = argparse.ArgumentParser()
    parser.add_argument("--codex", required=True, help="registered CODEX OME-TIFF")
    parser.add_argument("--he", default=None, help="registered H&E OME-TIFF")
    parser.add_argument("--markers", default=None, help="one channel name per line")
    parser.add_argument("--celltype-mask", default=None)
    parser.add_argument("--cn-mask", default=None)
    args = parser.parse_args()

    channel_names = (list(pd.read_csv(args.markers, header=None)[0])
                     if args.markers else None)
    labs = tifffile.imread(args.celltype_mask) if args.celltype_mask else None
    cn_labs = tifffile.imread(args.cn_mask) if args.cn_mask else None

    lazy_view(args.codex, channel_names=channel_names,
              labs=labs, cols=CELLTYPE_COLORS,
              cn_labs=cn_labs, cn_cols=CN_COLORS, he_path=args.he)
    napari.run()
