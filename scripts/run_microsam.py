import skimage
from datetime import datetime
import numpy as np
import tifffile as tf

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import warnings
from typing import Optional, Tuple

from skimage.filters import gaussian
from skimage.draw import disk
from skimage.io import imshow
#from scipy.ndimage import distance_transform_edt
from skimage.measure import regionprops
from skimage.morphology import binary_closing


import micro_sam.util as util
from micro_sam.prompt_based_segmentation import segment_from_points, segment_from_mask, segment_from_box, segment_from_box_and_points
from segment_anything.predictor import SamPredictor
from segment_anything.utils.transforms import ResizeLongestSide
import glob
import pickle
import skimage
from tqdm import tqdm
import os
#30 seconds
start = datetime.now(); print(start)

this_input = list(snakemake.input)
this_output = snakemake.output[0]

input_cell_body = this_input[0]
input_segmentation_for_seeds = this_input[1]

#print(f'{input_cell_body=}', f'{input_segmentation_for_seeds=}')
#print(f'{this_output=}')

labels = skimage.io.imread(input_segmentation_for_seeds)

centroid_dict = skimage.measure.regionprops_table(labels, properties=('centroid',))
#plt.imshow(labels)
#plt.scatter(centroid_dict['centroid-1'], centroid_dict['centroid-0'], s=2)

model_type = "vit_h"
predictor = util.get_sam_model(model_type=model_type)
print('model loaded: ', model_type)

def get_close_points(i, points, distance_threshold):
    """
    This function finds the points closest to the point indexed by i and returns them 
    """
    this_point = points[i]
    dist = np.linalg.norm(points-this_point, axis=1)
    indices_of_close_points = list(np.where(dist < distance_threshold)[0])
    indices_of_close_points.remove(i)
    return points[indices_of_close_points]

image = skimage.io.imread(input_cell_body)
image_embeddings = util.precompute_image_embeddings(predictor, image)
util.set_precomputed(predictor, image_embeddings)   

points = np.array([[x0, x1] for x1, x0 in zip(centroid_dict['centroid-1'],
                                   centroid_dict['centroid-0'])])
output_labels = np.zeros(image.shape, dtype=np.uint16)

j = 1
for i in range(len(points)):
    this_positive_point = points[i]
    these_negative_points = get_close_points(i, points, 0.0001)
    these_points = np.vstack((this_positive_point, these_negative_points))
    sam_labels = np.zeros(len(these_points))
    sam_labels[0] = 1

    predicted_points = segment_from_points(predictor, these_points, sam_labels)
    if len(output_labels[predicted_points[0] == 1])/output_labels.size < 0.5:
        output_labels[predicted_points[0] == 1] = j
        j += 1
    
    if False:
        fig, ax = plt.subplots(1,3, figsize=(13,4.5))
        ax[0].imshow(labels)
        ax[1].imshow(image)
        im = ax[2].imshow(output_labels)

        ax[0].scatter(this_positive_point[1], this_positive_point[0])
        ax[1].scatter(this_positive_point[1], this_positive_point[0])
        ax[2].scatter(this_positive_point[1], this_positive_point[0])
        
        ax[0].scatter(these_negative_points[:,1], these_negative_points[:,0])
        ax[1].scatter(these_negative_points[:,1], these_negative_points[:,0])
        ax[2].scatter(these_negative_points[:,1], these_negative_points[:,0])

        plt.colorbar(im, ax=ax[2])

        plt.show()
#    if i == 10:
#        break

skimage.io.imsave(this_output, output_labels, check_contrast=False)