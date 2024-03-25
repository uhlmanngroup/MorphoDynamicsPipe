import os, sys
import numpy as np
import skimage
import pandas as pd
import scipy
import natsort
from skimage.measure import label, regionprops, regionprops_table
import math

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


myprops = ('label', 'area', 'area_bbox', 'area_convex', 'area_filled',
            'axis_major_length', 'axis_major_length', 'eccentricity', 'equivalent_diameter_area',
           'euler_number', 'extent', 'feret_diameter_max', 'moments_hu', 'perimeter',
           'perimeter_crofton', 'solidity', 'bbox', 'centroid')

list_of_dataframes = []

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

def get_instantaneous_speeds(this_cell_df, cycleTime):
    # This function leaves the last value as a -1
    x_positions = this_cell_df['centroid-0'].values
    y_positions = this_cell_df['centroid-1'].values
    times = (this_cell_df.index*cycleTime).values
    speed = -1*np.ones(len(x_positions), float)
    speed_angle = -1*np.ones(len(x_positions), float)
    try:
        for i in range(len(x_positions) -1):
            displacement = (x_positions[i+1] - x_positions[i], 
                            y_positions[i+1] - y_positions[i])
            dist = np.linalg.norm(displacement)
            time_difference = times[i+1] - times[i]
            speed[i] = dist/time_difference
            possible_angle = math.degrees(np.arctan2(displacement[1], displacement[0]))
            # this part converts the angle to be between 0-360 degrees
            if possible_angle > 0:
                speed_angle[i] = possible_angle
            else:
                speed_angle[i] = 360 + possible_angle          
        return speed, speed_angle
    except:
        return speed, speed_angle
    
cell_indices = list(df_all.index.get_level_values(0).unique())

for this_cellID in cell_indices:
#    if this_cellID != 3:
#        continue
    this_cell_df = df_all.loc[this_cellID]
    speed, speed_angle = get_instantaneous_speeds(this_cell_df, cycleTime)
    df_all.loc[this_cellID, 'speed'] = speed
    df_all.loc[this_cellID, 'speed_angle'] = speed_angle
    df_all.loc[this_cellID, 'N_frames_existence'] = int(len(this_cell_df))
#    display(this_cell_df)

df_all.to_csv(this_output)