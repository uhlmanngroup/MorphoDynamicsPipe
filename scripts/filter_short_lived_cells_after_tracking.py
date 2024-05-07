#this file filters objects that do no last for a long time
import os
import skimage
import numpy as np
this_input = list(snakemake.input)
this_output = list(snakemake.output)

im = skimage.io.imread(this_input[0])
track_info= np.load(this_input[1])
tracked_ids = track_info[:,0].astype(int)
number_of_frames_to_remove_object = 8 #max number of consecutive frames for an object to be removed
ids_to_remove = np.bincount(tracked_ids) <= number_of_frames_to_remove_object
ids_in_current_image = np.unique(im)
for each_id_in_current_image in ids_in_current_image:
    if each_id_in_current_image in tracked_ids:
        if ids_to_remove[each_id_in_current_image]:
            im[im == each_id_in_current_image] = 0
        else:
            continue
    else:
        im[im == each_id_in_current_image] = 0
skimage.io.imsave(this_output[0], im, check_contrast=False)