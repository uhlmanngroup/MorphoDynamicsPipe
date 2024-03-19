import os, sys
import pickle
from glob import glob
import numpy as np
import skimage
from functools import partial
#import useful_functions as uf
import matplotlib.pyplot as plt
from tqdm import tqdm, trange

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
this_output = snakemake.output[0]

#print("this_input is in the python file for different inputs ", this_input)
#print('type ', type(this_input))
#print("this_output is in the python file for different inputs ", this_output)

frame = skimage.io.imread(this_input[0])
data = np.load(this_input[1])
T = getvaluefromstringbest(this_input[0], '_T=', preceding = '', ending='.tif', mydtype=int)
#print('T = ', T)

data_for_this_frame = data[data[:,1] == T]
dict_data_id_to_coords = {int(each[0]):(each[2], each[3]) for each in data_for_this_frame}
output_segmentation = np.zeros(frame.shape, frame.dtype)

for this_id in dict_data_id_to_coords.keys():
    this_coords = dict_data_id_to_coords[this_id]
    segmentation_id = frame[int(this_coords[0]), int(this_coords[1])]
    if segmentation_id == 0:
        continue
    output_segmentation += (frame == segmentation_id)*this_id

skimage.io.imsave(this_output, output_segmentation, check_contrast=False)