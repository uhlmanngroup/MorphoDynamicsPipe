import os, sys
import numpy as np
import pandas as pd
import skimage.io
from skimage.measure import label, regionprops, regionprops_table

this_input = sys.argv[1]
this_output = sys.argv[2]

def remove_labels_on_the_edge(masks):
    top = masks[0,:]
    bottom = masks[-1, :]
    left = masks[:, 0]
    right = masks[:, -1]
    myset = set()
    myset.update(top)
    myset.update(bottom)
    myset.update(left)
    myset.update(right)
    myset.remove(0)
    masks_out = masks.copy()
    for each_label in list(myset):
        masks_out[masks_out == each_label] = 0
    return masks_out


masks = skimage.io.imread(this_input)
masks_edge_removed0 =remove_labels_on_the_edge(masks)

myprops = ('label', 'area', 'area_bbox', 'area_convex', 'area_filled',
            'axis_major_length', 'axis_major_length', 'eccentricity', 'equivalent_diameter_area',
           'euler_number', 'extent', 'feret_diameter_max', 'moments_hu', 'perimeter',
           'perimeter_crofton', 'solidity')
#I removed num_pixels because it is not supported in some versions of skimage

# maybe inertia tensor is useful
df = pd.DataFrame(regionprops_table(masks_edge_removed0 , properties = myprops))

df.to_csv(this_output, header=True, index=False)

