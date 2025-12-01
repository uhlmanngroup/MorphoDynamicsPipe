from datetime import datetime
import numpy as np
import time, os, sys
from urllib.parse import urlparse
import skimage.io
from snakemake.script import snakemake

from urllib.parse import urlparse
from cellpose import models, core, denoise, utils
from cellpose.io import logger_setup

# logger_setup();

use_GPU = core.use_gpu()

# this_input = sys.argv[1]
# this_output = sys.argv[2]

this_input = list(snakemake.input)
print(this_input)
this_output = list(snakemake.output)

myparams = {
    "model_type": "cyto3",
    "restore_type": "denoise_cyto3",  #:None,
    "diameter": 30.0,
    "flow_threshold": 0.4,
    "cellprob_threshold": 0,
    "normalize": {"percentile": [1, 99]},
    "pretrained_model": "/example_file",
}

# help(snakemake)
# myparams = list(snakemake.params)[0]
# print(myparams)
# myparams = snakemake.Params
# print(snakemake.myparam0)
# print(snakemake.params.myparam0)
# print(snakemake.params['myparam0'])
# print(params)


im = skimage.io.imread(this_input[0])

# path_to_pretrained_model = os.path.abspath('../../2024-04-17_retrain_cellpose/
# 2D_nuclei_confocal_no_preprocess/OLD_test/models/CP_20240422_135130')

#<<<<<<< Updated upstream
#model = denoise.CellposeDenoiseModel(gpu=use_GPU, 
#                                     model_type=myparams['model_type'],
#                                     restore_type=myparams['restore_type'],
###                                     pretrained_model=myparams['pretrained_model'],
#)


#model = denoise.CellposeDenoiseModel(gpu=use_GPU, 

model = denoise.CellposeDenoiseModel(
    gpu=use_GPU,
    model_type=myparams["model_type"],
    restore_type=myparams["restore_type"],
    ##                                     pretrained_model=myparams['pretrained_model'],
)

labels, _, _, _ = model.eval(
    im,
    diameter=myparams["diameter"],
    flow_threshold=myparams["flow_threshold"],
    channels=[0, 0],
    cellprob_threshold=myparams["cellprob_threshold"],
    normalize=myparams["normalize"],
)

# labels, _, _, _ = model.eval(im, diameter=diameter, flow_threshold=flow_threshold, channels=[0, 0],
#    cellprob_threshold=cellprob_threshold, normalize=normalize)

skimage.io.imsave(this_output[0], labels, check_contrast=False)
# skimage.io.imsave(this_output[0], np.array([0]), check_contrast=False)
