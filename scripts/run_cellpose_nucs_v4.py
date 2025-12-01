from datetime import datetime
import numpy as np
import time, os, sys
from urllib.parse import urlparse
import skimage.io

# Patch torch before importing Cellpose to avoid BFloat16 on MPS
#import torch
#try:
#    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
#        # Map bfloat16 to float32 so .to(dtype) calls won't try to use bfloat16 on MPS
#        if hasattr(torch, "bfloat16"):
#            torch.bfloat16 = torch.float32
##        try:
#            torch.set_default_dtype(torch.float32)
#        except Exception:
#            pass
#except Exception:
#    pass

# now import cellpose after the torch patch
from cellpose import models, core, denoise, utils

# logger_setup();

use_GPU = core.use_gpu()

this_input = sys.argv[1]
this_output = sys.argv[2]

#this_input = list(snakemake.input)
print(this_input)
#this_output = list(snakemake.output)

my_params = {
    "model_type": "cyto3",
    "restore_type": "denoise_cyto3",  #:None,
    "diameter": 30.0,
    "flow_threshold": 0.4,
    "cellprob_threshold": 0,
    "normalize": {"percentile": [1, 99]},
    "pretrained_model": "cpsam",
    "max_iteration": 250,
}

im = skimage.io.imread(this_input)

# path_to_pretrained_model = os.path.abspath('../../2024-04-17_retrain_cellpose/
# 2D_nuclei_confocal_no_preprocess/OLD_test/models/CP_20240422_135130')

model = models.CellposeModel(
    gpu=use_GPU,
    pretrained_model=my_params["pretrained_model"],
    # model_type=myparams["model_type"],#deprecated
    # restore_type=myparams["restore_type"],#deprecated
    ##                                     pretrained_model=myparams['pretrained_model'],
)

labels, _, _= model.eval(
    im,
    niter=my_params["max_iteration"],
    diameter=my_params["diameter"],
    flow_threshold=my_params["flow_threshold"],
    # channels=[0, 0], #deprecated
    cellprob_threshold=my_params["cellprob_threshold"],
    # normalize=my_params["normalize"],#deprecated
)
# labels, _, _, _ = model.eval(im, diameter=diameter, flow_threshold=flow_threshold, channels=[0, 0],
#    cellprob_threshold=cellprob_threshold, normalize=normalize)
skimage.io.imsave(this_output, labels, check_contrast=False)
# skimage.io.imsave(this_output[0], np.array([0]), check_contrast=False)
