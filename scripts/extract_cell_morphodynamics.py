import os, sys
import numpy as np
import skimage
import pandas as pd
import scipy
import natsort

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
unique_ids = np.unique(data[:,0])

def speeds_weighted(this_cell_data, cycleTime):
    speeds = []
    for i in range(len(this_cell_data) - 1):
        displacement = this_cell_data[i+1, 2:4] - this_cell_data[i, 2:4]
        dist = np.linalg.norm(displacement)
        time_steps = this_cell_data[i+1, 1] - this_cell_data[i, 1]
        speeds += [dist/(time_steps*cycleTime)]*int(time_steps)
        mean_speed = np.mean(speeds)
        stdev_speed = np.std(speeds)
        min_speed = np.min(speeds)
        max_speed = np.max(speeds)
        median_speed = np.median(speeds)
        iqr_speed = scipy.stats.iqr(speeds)
    return mean_speed, stdev_speed, min_speed, max_speed, median_speed, iqr_speed

def time_span_properties(this_cell_data, cycleTime):
    start_frame = int(np.min(this_cell_data[:,1]))
    end_frame = int(np.max(this_cell_data[:,1]))
    N_separate_frames = len(this_cell_data)
    N_time_span_frames = end_frame - start_frame
    time_span = N_time_span_frames*cycleTime
    return start_frame, end_frame, N_separate_frames, N_time_span_frames, time_span

def get_area_from_segmentation(this_cell_seg):
    myareas = np.sum(this_cell_seg, axis=(1,2))
    myareas_non_zero = [each for each in myareas if each != 0]
    
    try:
        mean_area = np.mean(myareas_non_zero)
        stdev_area = np.std(myareas_non_zero)
        min_area = np.min(myareas_non_zero)
        max_area = np.max(myareas_non_zero)
        median_area = np.median(myareas_non_zero)
        iqr_area = scipy.stats.iqr(myareas_non_zero)

        return mean_area, stdev_area, min_area, max_area, median_area, iqr_area
    except:
        return -1, -1, -1, -1, -1, -1

list_cell_properties = []

for this_id in unique_ids:
    this_cell_data = data[data[:,0] == this_id]
    
    if len(this_cell_data) <= 1:
        continue
 #   print(this_cell_data)
    
    this_cell_properties = {'id':int(this_id)}
    
    start_frame, end_frame, N_separate_frames,\
    N_time_span_frames, time_span = time_span_properties(this_cell_data, cycleTime)
    this_cell_properties['start_frame'] = start_frame
    this_cell_properties['end_frame'] = end_frame
    this_cell_properties['N_separate_frames'] = N_separate_frames
    this_cell_properties['N_time_span_frames'] = N_time_span_frames
    this_cell_properties['time_span'] = time_span
    
    mean_speed, stdev_speed, min_speed,\
    max_speed, median_speed, iqr_speed = speeds_weighted(this_cell_data, cycleTime)
    this_cell_properties['mean_speed'] = mean_speed
    this_cell_properties['stdev_speed'] = stdev_speed
    this_cell_properties['min_speed'] = min_speed
    this_cell_properties['max_speed'] = max_speed
    this_cell_properties['median_speed'] = median_speed
    this_cell_properties['iqr_speed'] = iqr_speed
    
    this_cell_seg = seg_relabeled == this_id
    
    mean_area, stdev_area, min_area, \
    max_area, median_area, iqr_area = get_area_from_segmentation(this_cell_seg)
    this_cell_properties['mean_area'] = mean_area
    this_cell_properties['stdev_area'] = stdev_area
    this_cell_properties['min_area'] = min_area
    this_cell_properties['max_area'] = max_area
    this_cell_properties['median_area'] = median_area
    this_cell_properties['iqr_area'] = iqr_area
    

    list_cell_properties.append(this_cell_properties)
    
 #   if this_id == 4:
 #       break

#pd.DataFrame([1]).to_csv(this_output)
df = pd.DataFrame(list_cell_properties)
df.set_index('id', drop=True, inplace=True)
df.to_csv(this_output)