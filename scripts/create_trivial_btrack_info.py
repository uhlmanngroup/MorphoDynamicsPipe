import os, sys
import pickle
from glob import glob
import numpy as np
import skimage
from functools import partial
#import useful_functions as uf
import matplotlib.pyplot as plt
from tqdm import tqdm, trange

this_input = list(snakemake.input)
this_output = snakemake.output[0]
#this_input = [os.path.abspath('1_data/date_20250327_cellline_PDX67_condition_tissue_ebiid_70_C=0_ChName_DAPI_series_1PDX67_1-DAPI___seriesid_1_resolution_0pt1625_cycleTime_0/ebiid_70_C=0_seriesid_1_T=000.tif'),]
#this_output = os.path.abspath('testX.npy')
print(len(this_input),)
segmentation = skimage.io.imread(this_input[0])

myprops = skimage.measure.regionprops_table(segmentation, properties=['label', 'centroid', ])

output_data = np.zeros((len(myprops['label']), 4), dtype=int)
output_data[:, 0] = myprops['label']
output_data[:, 1] = 0  # Assuming all objects are in frame 0, as this is trivial
output_data[:, 2] = myprops['centroid-0']
output_data[:, 3] = myprops['centroid-1']
np.save(this_output, output_data)
