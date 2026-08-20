#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
from skimage import io
from matplotlib import pyplot as plt
from PIL import Image
import pandas as pd
import tifffile
import imageio
import glob
import re
import math
from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import make_outline_overlay
from deepcell.utils.plot_utils import create_rgb_image


# ---- Paths (edit these, or set the environment variables) -------------------

# CODEX_INPUT_DIR : directory holding the stacked tif, with a mesmer/ subfolder
# CODEX_SAMPLE    : sample identifier used in the stacked tif filename
inputDir = os.environ.get('CODEX_INPUT_DIR', 'CODEX/sample')
sampleId = os.environ.get('CODEX_SAMPLE', 'SAMPLE')

os.makedirs(os.path.join(inputDir, "mesmer"), exist_ok=True)

# Input OME.TIFF file and number of divisions for tiling
inOmeTiff = os.path.join(inputDir, sampleId + '_nuclear_panmembrane.tif')
ndivides = 3

# Reading and tiling image
print("Reading Image")
im = tifffile.imread(inOmeTiff)
print(im.shape)

print("Tiling Image")
M = im.shape[1] // ndivides
N = im.shape[2] // ndivides
tiles = [im[:, x:x + M, y:y + N] for x in range(0, im.shape[1], M) for y in range(0, im.shape[2], N)]

# Save tiled images
num = 1
for tile in tiles:
    suffix = str(num) + "_cell_panmembrane.qptiff"
    tifffile.imwrite(os.path.join(inputDir,"mesmer",suffix), tile)
    print(f"Printing tile: {num}")
    print(tile.shape)
    num += 1


# DeepCell access token -- set DEEPCELL_ACCESS_TOKEN in your environment.
# Request one at https://users.deepcell.org ; do not commit it to the repo.
if "DEEPCELL_ACCESS_TOKEN" not in os.environ:
    raise RuntimeError("DEEPCELL_ACCESS_TOKEN is not set in the environment")

# Function to run Mesmer on a tiled image, whole cell prediction
def run_mesmer(inputDir, filein, tilenum):
    outDir = inputDir + "/mesmer"

    print(f'Input OME.TIFF file: {filein}')
    outFileCell = outDir+"/"+ str(tilenum)+ "_segmentation_cell.tif"

    app = Mesmer()
    print('Training Resolution:', app.model_mpp, 'microns per pixel')
    im = io.imread(filein)
    X_train = np.transpose(im, (1, 2, 0))
    X_train = X_train.reshape(1, X_train.shape[0], X_train.shape[1], X_train.shape[2])
    rgb_images = create_rgb_image(X_train, channel_colors=['green', 'blue'])

    print('Starting Cell Segmentation')
    segmentation_predictions = app.predict(X_train, image_mpp=0.5)
    im = Image.fromarray(segmentation_predictions[0, :, :, 0]).save(outFileCell)
    overlay_data = make_outline_overlay(rgb_data=rgb_images, predictions=segmentation_predictions)
    fig, ax = plt.subplots(1, 2, figsize=(15, 15))
    ax[0].imshow(rgb_images[0, ...])
    ax[1].imshow(overlay_data[0, ...])
    ax[0].set_title('Raw data')
    ax[1].set_title('Dilated Cell Predictions')
    print("Saving overlay")
    fig.savefig(outDir+"/"+ str(tilenum)+ "_segmentation_cell_overlay.tif", dpi=600)
    with open(outDir+"/"+ str(tilenum)+ "_segmentation_predictions.npy", "wb") as outnp:
        np.save(outnp, segmentation_predictions)
    return

# Segmenting each tiled image using Mesmer
numimages = int(math.pow(ndivides, 2))
for j in range(1, numimages + 1):
    suffix = str(j) + "_cell_panmembrane.qptiff"
    imagetilepath = os.path.join(inputDir, "mesmer", suffix)
    print(f"Segmenting tile: {j}")
    print(os.listdir(inputDir+"/mesmer"))
    print(os.path.isfile(imagetilepath))
    run_mesmer(inputDir, imagetilepath, j)


# Merging whole cell numpy mask
print("Merging whole cell numpy mask")
extension = '_segmentation_predictions.npy'
imagetilearray = []
rowimagearray = []
curind = 0
for k in range(1, numimages + 1):
    print(curind)
    mask = np.load(inputDir +"/mesmer/"+str(k)+extension)
    maximum = mask.max()
    with np.nditer(mask, op_flags=['readwrite']) as it:
        for x in it:
            if x > 0:
                x[...] = x + curind
    imagetilearray.append(mask)
    curind = curind + maximum
ind = 0
for l in range(ndivides):
    x = ind
    y = ind + ndivides
    newrow = np.concatenate((imagetilearray[x:y]), axis=2)
    ind = ind + ndivides
    rowimagearray.append(newrow)
merged_whole_cell_mask = np.concatenate((rowimagearray), axis=1)
with open(inputDir+ "/mesmer/segmentation_predictions.npy", "wb") as outnp:
    np.save(outnp, merged_whole_cell_mask)

# Merging cell tiff mask
print("Merging whole cell tiff mask")
extension = '_segmentation_cell.tif'
imagetilearray = []
rowimagearray = []
curind = 0
for k in range(1, numimages + 1):
    print(curind)
    im = tifffile.imread(inputDir + "/mesmer/"+ str(k) + extension)
    maximum = im.max()
    with np.nditer(im, op_flags=['readwrite']) as it:
        for x in it:
            if x > 0:
                x[...] = x + curind
    imagetilearray.append(im)
    curind = curind + maximum
ind = 0
for l in range(ndivides):
    x = ind
    y = ind + ndivides
    newrow = np.concatenate((imagetilearray[x:y]), axis=1)
    ind = ind + ndivides
    rowimagearray.append(newrow)
merged_whole_cell_mask = np.concatenate((rowimagearray), axis=0)
tifffile.imwrite(inputDir + "/mesmer/segmentation_cell.tif", merged_whole_cell_mask)
