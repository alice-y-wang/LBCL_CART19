import numpy as np
from skimage import io
from matplotlib import pyplot as plt
#from deepcell.datasets import multiple tissue
#from deepcell.datasets import TissueNetSample
from deepcell.utils.plot_utils import create_rgb_image
from deepcell.utils.plot_utils import make_outline_overlay
from PIL import Image
import pandas as pd
import tifffile
import imageio
import glob
import multiprocessing as mp
import re
import os
Image.MAX_IMAGE_PIXELS = None

# ---- Paths (edit these, or set the environment variables) -------------------
# CODEX_DIR   : mcmicro output directory for one sample
# CODEX_SAMPLE: sample identifier used in the stacked tif filename
# CODEX_QPTIFF: raw multiplexed image
codexDir   = os.environ.get('CODEX_DIR', 'CODEX/sample')
sampleId   = os.environ.get('CODEX_SAMPLE', 'SAMPLE')
qptiffGlob = os.environ.get('CODEX_QPTIFF', 'CODEX/raw/*.qptiff')

mcmicroDir = codexDir
mcmicroDir2 = os.path.join(codexDir, 'mesmer')
inOmeTiff = glob.glob(qptiffGlob)[0]
print(inOmeTiff)

# Whole Cell
inMaskTiff = mcmicroDir2 + '/segmentation_cell.tif'
def sumValueByMask(mask, expr, outFile = ""):
    df = pd.DataFrame({
        'values':  np.concatenate(expr),
        'indexes': np.concatenate(mask)
    })
    out = df.pivot_table(values='values', index='indexes', aggfunc=sum)
    if(outFile != ""):
        out.to_csv(outFile)
    return out

mask = tifffile.imread(inMaskTiff)
print(mask.shape)
print(mask.dtype)
im = io.imread(os.path.join(mcmicroDir, sampleId + '_nuclear_panmembrane.tif'))
X_train = np.transpose(im, (1,2,0))
X_train = X_train.reshape(1, X_train.shape[0], X_train.shape[1], X_train.shape[2])
rgb_images = create_rgb_image(X_train, channel_colors=['green', 'blue'])
segmentation_predictions_nuc = np.load(mcmicroDir2 + "/segmentation_predictions.npy")
print(segmentation_predictions_nuc.shape)
mcmicroDir2 = os.path.join(codexDir, 'expr_matrix')
os.makedirs(mcmicroDir2, exist_ok=True)
df = pd.DataFrame({
  'values': np.concatenate(np.full(mask.shape,1)),
  'indexes': np.concatenate(mask)
})
size = sumValueByMask(mask, np.full(mask.shape,1), mcmicroDir2 + '/cell_size.csv')
coordX = sumValueByMask(
    mask,
    np.tile(
        np.array(range(0,mask.shape[1])),
        (mask.shape[0], 1)),
    mcmicroDir2 + '/cell_x_coord.csv')
coordY = sumValueByMask(
    mask,
    np.transpose(np.tile(
        np.array(range(0,mask.shape[0])),
        (mask.shape[1], 1))),
    mcmicroDir2 + '/cell_y_coord.csv')
im = tifffile.imread(inOmeTiff)
print(im.shape)
print(mask.shape)
print(im.max())
print(mask.max())
outFiles = [mcmicroDir2 + '/channel_C{:02d}.csv'.format(ii) for ii in range(im.shape[0])]
for ii in range(len(outFiles)):
    print(ii)
    sumValueByMask(mask, im[ii,:,:], outFiles[ii])

import os
csvDir = os.path.join(mcmicroDir2)
csvFileArray = sorted(glob.glob(csvDir + "/*.csv"))
print("Output Directory is " + csvDir)
df_list = []
for filename in csvFileArray:
    df = pd.read_csv(filename, header=0)
    df = df.rename(columns={"values": str(filename)}) # changes "values" to the filename
    df_list.append(df)
col_names = [str(i) for i in csvFileArray]
dfs = [df.set_index('indexes') for df in df_list]
frame = pd.concat(dfs, axis=1)
print(np.shape(frame)) # to confirm that the data-frame has written the correct number of columns
frame.to_csv(os.path.join(csvDir + '/combined_markers.csv'), index=True)
