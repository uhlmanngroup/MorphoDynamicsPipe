import numpy as np
import skimage
from numba import njit
from skimage.color import label2rgb

@njit
def get_edges_of_cluster_shapes_with_image_edges_int_output(intlabels, rgblabels, k = 2):
    """
    This function takes labels in two formats and returns the 
    label edges as a k-pixel edge image

    Parameters
    ----------
    intlabels : numpy array of ints with bg=0
        Labels in the format of an image with a different int for
        each label
    rgblabels : numpy array of floats RGB
        after from skimage.color import label2rgb has been called 
        on intlabels, with no image behind it
    k : int
        k is basically the thickness of the edging

    Returns
    -------
    output_edges 
        the image of output edges between the labels as int output

    """
    
    shape = intlabels.shape
    output_edges = np.zeros(shape, intlabels.dtype)
    for i in range(shape[0]):
        for j in range(shape[1]):
            this_val = intlabels[i,j]
            if this_val == 0:
                continue
            else:
                if i<k or shape[0]-i<k or j<k or shape[1]-j<k:
                    output_edges[i,j] = intlabels[i,j]
                else:
                    compare_array = intlabels[i-k:i+k+1,j-k:j+k+1]
                    allsame = np.all(compare_array == this_val)
                    
                    if allsame:
                        continue
                    else:
                        output_edges[i,j] = intlabels[i,j]
    return output_edges

this_input = list(snakemake.input)
this_output = list(snakemake.output)

labels = skimage.io.imread(this_input[0])
myedges = get_edges_of_cluster_shapes_with_image_edges_int_output(labels, label2rgb(labels), k=1)

skimage.io.imsave(this_output[0], myedges, check_contrast=False)
