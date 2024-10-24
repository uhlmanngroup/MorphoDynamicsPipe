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
myparams = list(snakemake.params)[0]

im = skimage.io.imread(this_input[0])

#path_to_pretrained_model = os.path.abspath('../../2024-04-17_retrain_cellpose/
#2D_nuclei_confocal_no_preprocess/OLD_test/models/CP_20240422_135130')

model = denoise.CellposeDenoiseModel(gpu=use_GPU, 
                                     model_type=myparams['model_type'],
             restore_type=myparams['restore_type'], 
#                                     pretrained_model=myparams['pretrained_model'],
                                    )

labels, _, _, _ = model.eval(im, diameter=myparams['diameter'], flow_threshold=myparams['flow_threshold'], channels=[0, 0],
    cellprob_threshold=myparams['cellprob_threshold'], normalize=myparams['normalize'])

skimage.io.imsave(this_output[0], labels, check_contrast=False)
