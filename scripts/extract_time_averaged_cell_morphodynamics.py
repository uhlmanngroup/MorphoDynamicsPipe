import os, sys
import numpy as np
import skimage
import pandas as pd
import scipy
import natsort

def getvaluefromstringbest(folder, variable, preceding='_', ending='_', mydtype=str):
    """
    This function is used to extract information from a string. 
    It is common used when a variable is encoded in the name of a file or folder.
    One example for this is 'Time_10_Channel_0.tif' 
    This function is used to extract the value of the variable 'Time' from the string with a value of 10

    Args:
    folder: string, this is the folder name (or any other string you would like to extract the value from)
    variable: string, this is the substring variable name that exists in the folder name
    preceding: string, this is the character between the variable and the value you would like to extract
    ending: string, this is the character that should follow the value you would like to extract
    If the value is at the end of the string, then I think it can be set to anything
    mydtype: type (e.g. int, str, float)
    """
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

data_name = this_input[0]
maximum_common_time_name = this_input[1]
segmentation_relabeled_names = this_input[2:]
image_names = [each.replace('3b_tracking_images', '1_data') for each in segmentation_relabeled_names]
try:
    cycleTime = getvaluefromstringbest(segmentation_relabeled_names[0], 
                                    'cycleTime', ending='/', mydtype=int)
except:
    cycleTime = 1

#print("tracking_info", tracking_info)
#print("segmentation", segmentation)

#print("this_output is in the python file for different inputs ", this_output)

data = pd.read_csv(data_name, index_col=(0))
data.sort_values(['cellID', 'frame_id_T'], inplace=True)

f = open(maximum_common_time_name, "r")
maximum_common_time = int(f.read())
f.close()

#seg_relabeled = np.array([skimage.io.imread(each) for each in seg_relabeled_names])
#images = np.array([skimage.io.imread(each) for each in image_names])

#This line limits the input data so that it works over a common time range
df_time = data.loc[data['realTime'] <= maximum_common_time]

def speeds_weighted(this_cell_data, cycleTime):
    speed_summaries = {}
    speeds = []
    myfuncs = [np.mean, np.std, np.min, np.max, np.median, scipy.stats.iqr]
    this_cell_data_np = this_cell_data[['frame_id_T', 'centroid-0', 'centroid-1']].to_numpy()
    try:
        for i in range(len(this_cell_data_np) - 1):
            displacement = this_cell_data_np[i+1, 1:3] - this_cell_data_np[i, 1:3]
            dist = np.linalg.norm(displacement)
            time_steps = this_cell_data_np[i+1, 0] - this_cell_data_np[i, 0]
            speeds += [dist/(time_steps*cycleTime)]*int(time_steps)
            for each_func in myfuncs:
                speed_summaries['speed2_' + each_func.__name__] = each_func(speeds)
        if len(speed_summaries.keys()) == 0:
            for each_func in myfuncs:
                speed_summaries['speed2_' + each_func.__name__] = -1
    except:
        for each_func in myfuncs:
            speed_summaries['speed2_' + each_func.__name__] = -1
    return speed_summaries

def time_span_properties(this_cell_data, cycleTime):
    dict_time_span_properties = {}
    if np.any(this_cell_data['frame_id_T'] == -1):
        dict_time_span_properties['start_frame'] = -1
        dict_time_span_properties['end_frame'] = -1
        dict_time_span_properties['N_separate_frames'] = -1
        dict_time_span_properties['N_time_span_frames'] = -1
        dict_time_span_properties['time_span'] = -1
    else:
        try:
            dict_time_span_properties['start_frame'] = int(np.min(this_cell_data['frame_id_T']))
            dict_time_span_properties['end_frame'] = int(np.max(this_cell_data['frame_id_T']))
            dict_time_span_properties['N_separate_frames'] = len(this_cell_data)
            dict_time_span_properties['N_time_span_frames'] = \
            dict_time_span_properties['end_frame'] - dict_time_span_properties['start_frame']
            dict_time_span_properties['N_missing_frames'] = \
            dict_time_span_properties['N_time_span_frames'] - dict_time_span_properties['N_separate_frames'] + 1 
            dict_time_span_properties['time_span'] = \
            dict_time_span_properties['N_time_span_frames'] * cycleTime
        except:
            dict_time_span_properties['start_frame'] = -1
            dict_time_span_properties['end_frame'] = -1
            dict_time_span_properties['N_separate_frames'] = -1
            dict_time_span_properties['N_time_span_frames'] = -1
            dict_time_span_properties['N_missing_frames'] = -1
            dict_time_span_properties['time_span'] = -1            
            
    return dict_time_span_properties

def get_area_from_segmentation(this_cell_seg):
    '''This just produces the same as regionprops so is not included
    '''
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
    
def get_summary_statistics(this_cell_data):
    columns = list(this_cell_data.columns)
    for each_remove in ['frame_id_T', 'realTime', 'frame_shape0', 'frame_shape1',
                        'speed', 'speed_angle', 'N_frames_existence']:
        columns.remove(each_remove)
    summary_stats = {}
    for each_col in columns:
        for function in [np.mean, np.std, np.min, np.max, np.median, scipy.stats.iqr]:
            try:
                if np.any(this_cell_data[each_col] == -1):
                    summary_stats[each_col + '_' + function.__name__] = -1
                else:
                    summary_stats[each_col + '_' + function.__name__] = function(this_cell_data[each_col])
            except:
                summary_stats[each_col + '_' + function.__name__] = -1
    return summary_stats


unique_cell_ids = np.unique(df_time.index.get_level_values(0))
list_cell_properties = []

#this interates over all the cells, calculating the morphodynamic properties
for this_cell_id in unique_cell_ids:
    this_cell_data= df_time.loc[df_time.index == this_cell_id]
    
    this_cell_properties = {'cellId':int(this_cell_id)}
    this_cell_properties.update(time_span_properties(this_cell_data, cycleTime))
    this_cell_properties.update(get_summary_statistics(this_cell_data))
    this_cell_properties.update(speeds_weighted(this_cell_data, cycleTime))
    
    #this is the XYT mask of the cell in the segmentation images
#    this_cell_seg = seg_relabeled == int(this_cell_id)
#    this_cell_images = images*this_cell_seg

    list_cell_properties.append(this_cell_properties)
    
#    if this_cell_id == 4:
#        break


df = pd.DataFrame(list_cell_properties)
df.set_index('cellId', drop=True, inplace=True)
df.to_csv(this_output)