import numpy as np
import os, sys
import skimage.io

from cellpose import models, core, denoise, utils
#from cellpose.io import logger_setup
#logger_setup();

use_GPU = core.use_gpu()

#this_input = sys.argv[1]
#this_output = sys.argv[2]

this_input = list(snakemake.input)
this_output = list(snakemake.output)
#myparams = list(snakemake.params)[0]

myparams = {'model_type':'cyto3',
    'restore_type':"denoise_cyto3",   #:None, 
    'diameter':30.0, 
    'flow_threshold':0.4, 
    'cellprob_threshold':0, 
    'normalize':{'percentile':[1, 99]},
    'pretrained_model':'/example_file',
    }

im_celltracker = skimage.io.imread(this_input[0])
im_nucs = skimage.io.imread(this_input[1])
im = np.array([im_celltracker, im_nucs, im_nucs])

#path_to_pretrained_model = os.path.abspath('../../2024-04-17_retrain_cellpose/
#2D_nuclei_confocal_no_preprocess/OLD_test/models/CP_20240422_135130')

model = denoise.CellposeDenoiseModel(gpu=use_GPU, 
                                     model_type=myparams['model_type'],
             restore_type=myparams['restore_type'], 
#                                     pretrained_model=myparams['pretrained_model'],
                                    )

labels, _, _, _ = model.eval(im, diameter=myparams['diameter'], flow_threshold=myparams['flow_threshold'], channels=[1, 2],
    cellprob_threshold=myparams['cellprob_threshold'], normalize=myparams['normalize'])
#Note, if you set the first channel input to use grayscale 0, 
#then no nuclear channel will be used (the second channel will be filled with zeros).

skimage.io.imsave(this_output[0], labels, check_contrast=False)