#uses view3 environment (no numba and napari=0.4.18)
import napari, os, skimage, pathlib
import numpy as np
from skimage.color import label2rgb
from magicgui import magicgui
from napari.types import LabelsData, ImageData, LayerDataTuple
from napari.layers import Layer, Labels, Image
from typing import List
import natsort
import dask_image.imread
import warnings
warnings.filterwarnings("ignore")

@magicgui(
    folder={"label": "Folder to open:\n1data/subfolder/", "mode": "d"},
    open_1_data={"widget_type": "CheckBox", "label": "1_data"},
    open_2_segmentation={"widget_type": "CheckBox", "label": "2_segmentation"},
    open_3b_tracking_images = {"widget_type": "CheckBox", "label": "3b_tracking_images"},
    open_3c_tracking_images_filtered = {"widget_type": "CheckBox", "label": "3c_tracking_images_filtered"},
    open_5_tracking_images_outlines = {"widget_type": "CheckBox", "label": "5_tracking_images_outlines"},
    lazy_loading={"widget_type": "CheckBox", "label": "Lazy loading"},
)
def open_folder(folder: pathlib.Path = '/nfs/research/uhlmann/', 
#                folder: pathlib.Path = '/Users/bwoodhams/Documents/plast_cell/1_data/',
                open_1_data:bool = True,
                open_2_segmentation:bool = False,
                open_3b_tracking_images:bool = False,
                open_3c_tracking_images_filtered:bool = False,
                open_5_tracking_images_outlines:bool = False,
                lazy_loading:bool = True,
                ) -> LayerDataTuple: #if put List here it doesn't work
    list_of_files = [each for each in os.listdir(folder) if each.endswith('.tif')]
    list_of_files = natsort.natsorted(list_of_files)

    dict_image_or_labels = {'1_data':'image', 
                            '2_segmentation':'labels', 
                            '3b_tracking_images':'labels', 
                            '3c_tracking_images_filtered':'labels', 
                            '5_tracking_images_outlines':'labels'}
    output_list = []

    if lazy_loading:
        for open_this, name in zip([open_1_data, open_2_segmentation, open_3b_tracking_images, open_3c_tracking_images_filtered, open_5_tracking_images_outlines], 
                                ['1_data', '2_segmentation', '3b_tracking_images', '3c_tracking_images_filtered', '5_tracking_images_outlines']):
            if open_this:
                this_folder = str(folder).replace('1_data', name)
                mylistofimages = dask_image.imread.imread(os.path.join(this_folder, '*.tif'))
                output_list.append((mylistofimages, {'name': name}, dict_image_or_labels[name]))
    else:
        for open_this, name in zip([open_1_data, open_2_segmentation, open_3b_tracking_images, open_3c_tracking_images_filtered, open_5_tracking_images_outlines], 
                                ['1_data', '2_segmentation', '3b_tracking_images', '3c_tracking_images_filtered', '5_tracking_images_outlines']):
            if open_this:
                this_folder = str(folder).replace('1_data', name)
                mylistofimages = np.array([skimage.io.imread(os.path.join(this_folder,file)) for file in list_of_files])
                output_list.append((mylistofimages, {'name': name}, dict_image_or_labels[name]))

    return output_list

viewer = napari.Viewer()
viewer.window.add_dock_widget(open_folder, area='bottom')
napari.run()
