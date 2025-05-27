#this file filters objects that do no last for a long time
import os, re
import skimage
import numpy as np
this_input = list(snakemake.input)
this_output = list(snakemake.output)
#number_of_frames_to_remove_object = snakemake.params[0] #max number of consecutive frames for an object to be removed
number_of_frames_to_remove_object = 0

im = skimage.io.imread(this_input[0])
track_info= np.load(this_input[1])
tracked_ids = track_info[:,0].astype(int)

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

# this part removes the entire object in the single time frame that at all overlaps with a correction labels
# file handling not currently included in the snakemake because it should be a dependency
correction_labels_name = re.sub('3b_tracking_images', '3b2_corrections', str(os.path.split(this_input[0])[0]))
correction_labels_name = os.path.join(correction_labels_name, 'correction_labels.tif')
print(correction_labels_name)
if os.path.exists(correction_labels_name):
    print('Activated')
    correction_labels = skimage.io.imread(correction_labels_name)
    correction_at_any_time = correction_labels == 1
    uniques = np.unique(im[correction_at_any_time])
    for each_unique in uniques:
        im[im == each_unique] = 0

# Save the filtered image
skimage.io.imsave(this_output[0], im, check_contrast=False)