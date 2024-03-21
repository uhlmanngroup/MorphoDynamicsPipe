import os, sys
import numpy as np
import skimage
import pandas as pd
import scipy
import natsort
from skimage.measure import label, regionprops, regionprops_table

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

tracking_info = this_input[0]
segmentation_relabeled_names = this_input[1:]
cycleTime = getvaluefromstringbest(segmentation_relabeled_names[0], 
                                   'cycleTime', ending='/', mydtype=int)

#print("tracking_info", tracking_info)
#print("segmentation", segmentation)

#print("this_output is in the python file for different inputs ", this_output)

data = np.load(tracking_info)
seg_relabeled = np.array([skimage.io.imread(each) for each in segmentation_relabeled_names])
unique_cell_ids = np.unique(data[:,0])
unique_frame_ids = np.unique(data[:,1])

list_of_dataframes = []

myprops = ('label', 'area', 'area_bbox', 'area_convex', 'area_filled',
            'axis_major_length', 'axis_major_length', 'eccentricity', 'equivalent_diameter_area',
           'euler_number', 'extent', 'feret_diameter_max', 'moments_hu', 'perimeter',
           'perimeter_crofton', 'solidity', 'bbox', 'centroid')

for this_frame_id in unique_frame_ids:
    this_frame = seg_relabeled[int(this_frame_id)]
    df = pd.DataFrame(regionprops_table(this_frame, properties = myprops))
    list_of_cellIDs = list(df['label'])
    df['frame_id_T'] = int(this_frame_id)
    df['realTime'] = int(this_frame_id) * cycleTime
    df['frame_shape0'] = this_frame.shape[0]
    df['frame_shape1'] = this_frame.shape[1]
    df2 = df.rename({'label':'cellID'}, axis=1).set_index(['cellID', 'frame_id_T'])
#    mycellid_props = df2.to_dict('index')
    
    this_frame_data = data[data[:,1] == this_frame_id]
    for this_cellID in list_of_cellIDs:
        this_cellID_frame_data = this_frame_data[this_frame_data[:,0] == this_cellID]
        if len(this_cellID_frame_data) == 1:
            df2.loc[this_cellID, 'track_pos0'] = this_cellID_frame_data[0, 2]
            df2.loc[this_cellID, 'track_pos1'] = this_cellID_frame_data[0, 3]
        else:
            df2.loc[this_cellID, 'track_pos0'] = -1
            df2.loc[this_cellID, 'track_pos1'] = -1
            
    list_of_dataframes.append(df2)

    df_all = pd.concat(list_of_dataframes)
    df_all.to_csv(this_output)