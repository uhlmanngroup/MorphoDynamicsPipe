import os, sys
import numpy as np
import btrack
import os
import pandas as pd
#from aicsimageio import AICSImage
import matplotlib.pyplot as plt
from tqdm import trange
import skimage
import numpy as np
#from skimage.io import imread
from glob import glob
#import useful_functions as uf
from functools import partial
#from btrack import datasets
import pickle

#this_input = sys.argv[1]
#this_output = sys.argv[2]

this_input = list(snakemake.input)
this_output = snakemake.output[0]
#myparams = list(snakemake.params)[0]
myparams = {'properties': (), #'area',
         'config_file': 'btrack_cell_config.json',
         'optimize': True,
        }

print('This_input', this_input)

#print("this_input is in the python file ", this_input)
#print('type ', type(this_input))

#print("this_output is in the python file ", this_output)

segmentation = np.array([skimage.io.imread(each) for each in this_input])
#print(segmentation.shape)

objects = btrack.utils.segmentation_to_objects(
  segmentation, properties=myparams['properties'],
)

with btrack.BayesianTracker() as tracker:

  # configure the tracker using a config file
    tracker.configure(myparams['config_file'])

  # append the objects to be tracked
    tracker.append(objects)

  # set the volume (Z axis volume limits default to [-1e5, 1e5] for 2D data)
    tracker.volume = ((0, segmentation.shape[2]), (0, segmentation.shape[1]))

  # track them (in interactive mode)
#  tracker.track_interactive(step_size=100)
    tracker.track(step_size = 100)
    # The number of tracking steps to be taken before returning summary
    # statistics. The tracking will be followed to completion, regardless
    # of the step size provided.

  # generate hypotheses and run the global optimizer
    if myparams['optimize']:
      tracker.optimize()
#    tracker.optimize()

  # store the data in an HDF5 file
 #   tracker.export(os.path.join(experiment_folder, 'hdf5_data.h5'), obj_type='obj_type_1')

  # get the tracks as a python list
    tracks = tracker.tracks

  # optional: get the data in a format for napari
    data, properties, graph = tracker.to_napari()


with open(this_output.replace('track_info.npy', 'properties.pkl'), "wb") as f:
    pickle.dump(properties, f)

with open(this_output.replace('track_info.npy', 'graph.pkl'), "wb") as f:
    pickle.dump(graph, f)

np.save(this_output, data)