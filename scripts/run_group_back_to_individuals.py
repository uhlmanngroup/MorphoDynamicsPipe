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


myarray = np.load(this_input[0], mmap_mode='r')
T = getvaluefromstringbest(this_input[1], '_T=', preceding='', ending='.tif', mydtype=int)


skimage.io.imsave(this_output[0], myarray[T], check_contrast=False)