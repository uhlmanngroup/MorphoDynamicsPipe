import os, sys
import numpy as np
import os
import pandas as pd
from tqdm import trange
import skimage
from functools import partial
from pystackreg import StackReg
import pickle as pkl

def getvaluefromstringbest(folder, variable, preceding='_', ending='_', mydtype=str):
    i = folder.index(variable)
    length = len(variable)
    i_ans_start = i + length + len(preceding)
    
    if ending in folder[i_ans_start:]:
        i_ans_end = folder[i_ans_start:].index(ending) + i_ans_start
    else:
        i_ans_end = len(folder)
    return mydtype(folder[i_ans_start:i_ans_end])

def normalize_frame_to_01(frame):
    mn, mx = frame.min(), frame.max()
    if mx > mn:
        return (frame - mn) / (mx - mn)
    return frame  # constant image edge-case

#this_input = sys.argv[1]
#this_output = sys.argv[2]

this_input = list(snakemake.input)
this_output = list(snakemake.output)

function_file = this_input[-1]
this_input = this_input[:-1]

#print("this_input is in the python file ", this_input)
#print('type ', type(this_input))

#print("this_output is in the python file ", this_output)

get_T = partial(getvaluefromstringbest, variable='_T=', preceding = '', ending='.tif', mydtype=int)
this_input.sort(key = get_T)
#this_output.sort(key = get_T)

image_stack = np.array([skimage.io.imread(each) for each in this_input])
#print(image_stack.shape)

if False: #normalizing to 0-1 before registration seems to help with the registration quality, but it is not strictly necessary
    image_stack = np.array([normalize_frame_to_01(frame) for frame in image_stack])
    print(image_stack.shape, image_stack.dtype)

#sr = StackReg(StackReg.TRANSLATION)
with open(function_file, 'rb') as f:
    sr = pkl.load(f)
out_previous = sr.transform_stack(image_stack)
#print(out_previous.shape)
#area_to_keep = np.all(out_previous != 0, axis=0)
#print(area_to_keep.shape)
#mytable = skimage.measure.regionprops_table(area_to_keep*1)
#print(mytable)
#cropped = out_previous[:, mytable['bbox-0'][0]:mytable['bbox-2'][0],
#             mytable['bbox-1'][0]:mytable['bbox-3'][0]]
#if image_stack.dtype == np.uint8:
#    cropped[cropped < 0] = 0
#    cropped[cropped > 255] = 255
#    cropped = cropped.astype(image_stack.dtype)

if image_stack.dtype == np.uint8:
    out_previous[out_previous < 0] = 0
    out_previous[out_previous > 255] = 255
    out_previous = out_previous.astype(image_stack.dtype)
if image_stack.dtype == np.uint16:
    out_previous[out_previous < 0] = 0
    out_previous[out_previous > 65535] = 65535
    out_previous = out_previous.astype(image_stack.dtype)

np.save(this_output[0], out_previous)
#np.save(this_output[0], cropped)

