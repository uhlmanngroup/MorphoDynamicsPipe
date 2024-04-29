from datetime import datetime
import numpy as np
import time, os, sys
from urllib.parse import urlparse
import skimage.io

from urllib.parse import urlparse
from cellpose import models, core, denoise, utils
from cellpose.io import logger_setup
#logger_setup();

use_GPU = core.use_gpu()

#this_input = sys.argv[1]
#this_output = sys.argv[2]

this_input = list(snakemake.input)
this_output = list(snakemake.output)

im = skimage.io.imread(this_input[0])

#path_to_pretrained_model = os.path.abspath('../../2024-04-17_retrain_cellpose/
#2D_nuclei_confocal_no_preprocess/OLD_test/models/CP_20240422_135130')

model = denoise.CellposeDenoiseModel(gpu=use_GPU, 
                                     model_type="cyto3",
             restore_type="denoise_cyto3", 
#                                     pretrained_model=path_to_pretrained_model,
                                    )

labels, _, _, _ = model.eval(im, diameter=20.0, flow_threshold=0.7, channels=[0, 0],
    cellprob_threshold=-2, normalize={'percentile':[0, 100]})

skimage.io.imsave(this_output[0], labels, check_contrast=False)
