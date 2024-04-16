import os, sys
import numpy as np
import os
import pandas as pd
from tqdm import trange
import skimage
from functools import partial
from pystackreg import StackReg
from datetime import datetime

def getvaluefromstringbest(folder, variable, preceding='_', ending='_', mydtype=str):
    i = folder.index(variable)
    length = len(variable)
    i_ans_start = i + length + len(preceding)
    
    if ending in folder[i_ans_start:]:
        i_ans_end = folder[i_ans_start:].index(ending) + i_ans_start
    else:
        i_ans_end = len(folder)
    return mydtype(folder[i_ans_start:i_ans_end])

#this_input = sys.argv[1]
#this_output = sys.argv[2]

this_input = list(snakemake.input)
this_output = list(snakemake.output)

print("this_input is in the python file ", this_input)
print('type ', type(this_input))

print("this_output is in the python file ", this_output)

get_T = partial(getvaluefromstringbest, variable='_T=', preceding = '', ending='.tif', mydtype=int)
this_input.sort(key = get_T)
#this_output.sort(key = get_T)
image_stack = np.array([skimage.io.imread(each) for each in this_input])

sr = StackReg(StackReg.TRANSLATION)
out_previous = sr.register_transform_stack(image_stack, reference='previous')
out_previous = out_previous.astype(image_stack.dtype)
area_to_keep = np.all(out_previous != 0, axis=0)
mytable = skimage.measure.regionprops_table(area_to_keep*1)
cropped = out_previous[:, mytable['bbox-0'][0]:mytable['bbox-2'][0],
             mytable['bbox-1'][0]:mytable['bbox-3'][0]]

#for i, this_output in enumerate(this_output):
#    skimage.io.imsave(this_output, cropped[i], check_contrast=False)
#    skimage.io.imsave(this_output, np.array([0]), check_contrast=False)
#skimage.io.imsave(this_output[0], cropped[i], check_contrast=False)
np.save(this_output[0], cropped)