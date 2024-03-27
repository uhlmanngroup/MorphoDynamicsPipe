#this file takes a list of folders as input, extracts their subfiles
# and finds the maximum time that is common to all time series
#This is to enable

import os
import natsort

this_input = list(snakemake.input)
this_output = snakemake.output[0]

def getvaluefromstringbest(folder, variable, preceding='_', ending='_', mydtype=str):
    i = folder.index(variable)
    length = len(variable)
    i_ans_start = i + length + len(preceding)
    
    if ending in folder[i_ans_start:]:
        i_ans_end = folder[i_ans_start:].index(ending) + i_ans_start
    else:
        i_ans_end = len(folder)
    return mydtype(folder[i_ans_start:i_ans_end])


realTime_common = 1000000000000000000000000000000

for each_subfolder in this_input:
    list_of_files = natsort.natsorted(os.listdir(each_subfolder))
    maxT = getvaluefromstringbest(list_of_files[-1], '_T=', preceding='', 
                           ending='.tif', mydtype=int)
    cycleTime = getvaluefromstringbest(each_subfolder, 'cycleTime', preceding='_', 
                           mydtype=int)
    realTime_max = maxT * cycleTime
    if realTime_max < realTime_common:
        realTime_common = realTime_max

f = open(this_output, "w")
f.write(str(realTime_common))
f.close()