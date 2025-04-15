#uses view3 environment (no numba and napari=0.4.18)
#rerun from 3c onwards after this 
import napari, os, skimage, pathlib, re
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
    folder2={"label": "Folder2 to open:\ncelltracker/1_data/subfolder/", "mode": "d"},
    folder3={"label": "Folder3 to open:\nbrightfield/1_data/subfolder/", "mode": "d"},
    open_1_data={"widget_type": "CheckBox", "label": "1_data"},
    open_3b_tracking_images = {"widget_type": "CheckBox", "label": "3b_tracking_images"},
    open_folder2 = {"widget_type": "CheckBox", "label": "Open folder2"},
    open_folder3 = {"widget_type": "CheckBox", "label": "Open folder3"},
    lazy_loading={"widget_type": "CheckBox", "label": "Lazy loading"},
)
def open_folder(folder: pathlib.Path = '/Users/bwoodhams/nfsresearch/bwoodhams/plast_cell/2024-08-22_making_example/MorphoDynamicsPipe/1_data/fieldofview1', 
#                folder: pathlib.Path = '/nfs/research/uhlmann/',
                folder2: pathlib.Path = '/Users/bwoodhams/nfsresearch/bwoodhams/plast_cell/2024-08-22_making_example/MorphoDynamicsPipe/1_data/fieldofview1', 
                folder3: pathlib.Path = '/Users/bwoodhams/nfsresearch/bwoodhams/plast_cell/2024-08-22_making_example/MorphoDynamicsPipe/1_data/fieldofview1', 
                open_1_data:bool = True,
                open_3b_tracking_images:bool = False,
                open_folder2:bool = False,
                open_folder3:bool = False,
                lazy_loading:bool = True,
                ) -> LayerDataTuple: #if put List here it doesn't work
    list_of_files = [each for each in os.listdir(folder) if each.endswith('.tif')]
    list_of_files = natsort.natsorted(list_of_files)

    dict_image_or_labels = {'1_data':'image', 
                            '3b_tracking_images':'labels', 
}
    output_list = []

    if lazy_loading:
        for open_this, name in zip([open_1_data, open_3b_tracking_images, ], 
                                ['1_data', '3b_tracking_images',]):
            if open_this:
                this_folder = str(folder).replace('1_data', name)
                mylistofimages = dask_image.imread.imread(os.path.join(this_folder, '*.tif'))
                output_list.append((mylistofimages, {'name': name}, dict_image_or_labels[name]))
        if open_folder2:
            this_folder = folder2
            mylistofimages2 = dask_image.imread.imread(os.path.join(this_folder, '*.tif'))
            output_list.insert(1, (mylistofimages2, {'name': 'folder2', 
                                                     'blending': 'additive',
                                                     'colormap': 'green'}, 'image', ))
        if open_folder3:
            this_folder = folder3
            mylistofimages3 = dask_image.imread.imread(os.path.join(this_folder, '*.tif'))
            output_list.insert(0, (mylistofimages3, {'name': 'folder3', 
                                                     'blending': 'additive',
                                                     'colormap': 'gray_r'}, 'image', ))
    else:
        for open_this, name in zip([open_1_data, open_3b_tracking_images, ], 
                                ['1_data', '3b_tracking_images', ]):
            if open_this:
                this_folder = str(folder).replace('1_data', name)
                mylistofimages = np.array([skimage.io.imread(os.path.join(this_folder,file)) for file in list_of_files])
                output_list.append((mylistofimages, {'name': name}, dict_image_or_labels[name]))
        if open_folder2:
            this_folder = folder2
            list_of_files2 = [each for each in os.listdir(this_folder) if each.endswith('.tif')]
            list_of_files2 = natsort.natsorted(list_of_files2)
            mylistofimages2 = np.array([skimage.io.imread(os.path.join(this_folder,file)) 
                                        for file in list_of_files2])
            output_list.insert(1, (mylistofimages2, {'name': 'folder2', 
                                                     'blending':'additive',
                                                     'colormap': 'green',
                                                     }, 'image',))
        if open_folder3:
            this_folder = folder3
            list_of_files3 = [each for each in os.listdir(this_folder) if each.endswith('.tif')]
            list_of_files3 = natsort.natsorted(list_of_files3)
            mylistofimages3 = np.array([skimage.io.imread(os.path.join(this_folder,file)) 
                                        for file in list_of_files3])
            output_list.insert(0, (mylistofimages3, {'name': 'folder3', 
                                                     'blending':'additive',
                                                     'colormap': 'gray_r'}, 'image',))

    # Add the correction labels layer
    correction_labels_name = re.sub('1_data', '3b2_corrections', str(folder))
    correction_labels_name = os.path.join(correction_labels_name, 'correction_labels.tif')
    if os.path.exists(correction_labels_name):
        correction_labels = skimage.io.imread(correction_labels_name)
    else:
        correction_labels = np.zeros(mylistofimages[0].shape, dtype=np.uint8)
    output_list.append((correction_labels, {'name': 'correction_labels'}, 'labels',))
    return output_list

@magicgui(To={"widget_type": "Label", "value": "save, click Run"})
def save_func(To, viewer: napari.Viewer):
    # Save the correction labels
    correction_labels = viewer.layers['correction_labels'].data
    data_folder = open_folder.folder.value
    output_folder = re.sub('1_data', '3b2_corrections', str(data_folder))
    os.makedirs(output_folder, exist_ok=True)
    skimage.io.imsave(os.path.join(output_folder, 'correction_labels.tif'), 
                      correction_labels,
                      check_contrast=False,)

viewer = napari.Viewer()
viewer.window.add_dock_widget(open_folder, area='bottom')
viewer.window.add_dock_widget(save_func, area='bottom')
napari.run()
