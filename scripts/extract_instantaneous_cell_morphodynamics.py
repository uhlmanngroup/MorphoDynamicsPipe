import os, sys
import numpy as np
import skimage
import pandas as pd
import scipy
import natsort
from skimage.measure import label, regionprops, regionprops_table
import math
import scipy.ndimage as nd
# print(dir(scripts))
# import scripts.cpda as cpda
import cpda as cpda
from functools import partial
from cp_measure.bulk import get_core_measurements

def get_number_of_neighbors(i, labels, pixel_range):
    mask = labels == i
    dilated = mask.copy()
    for i in range(pixel_range):
        dilated = skimage.morphology.binary_dilation(dilated)
    edge = dilated ^ mask
    neighbors = np.unique(labels[edge])
    if 0 in neighbors:
        return len(neighbors) - 1
    else:
        return len(neighbors)

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

tracking_info = this_input.pop(0)
segmentation_relabeled_names = this_input[0:int(len(this_input)/2)]
images_names = this_input[int(len(this_input)/2):]

try:
    cycleTime = getvaluefromstringbest(segmentation_relabeled_names[0], 
                                    'cycleTime', ending='/', mydtype=int)
except:
    cycleTime = 1

#print("tracking_info", tracking_info)
#print("segmentation", segmentation)

#print("this_output is in the python file for different inputs ", this_output)

data = np.load(tracking_info)
seg_relabeled = np.array([skimage.io.imread(each) for each in segmentation_relabeled_names])
images = np.array([skimage.io.imread(each) for each in images_names])

if False:
    #this part should be set to true for testing on a small part of the dataset
    #otherwise should be false
    list_of_images_to_keep = [0, 1, 2, 3, 4]
    myfilter = np.zeros(len(data[:,1]), np.bool_)
    for i in [0, 1, 2, 3, 4]:
        myfilter = np.logical_or(data[:,1] == i, myfilter)
    data = data[myfilter]
    seg_relabeled = seg_relabeled[list_of_images_to_keep]
    images = images[list_of_images_to_keep]

unique_cell_ids = np.unique(data[:,0])
unique_frame_ids = np.unique(data[:,1])


myprops = ('label', 'area', 'area_bbox', 'area_convex', 'area_filled',
            'axis_major_length', 'axis_major_length', 'eccentricity', 'equivalent_diameter_area',
           'euler_number', 'extent', 'feret_diameter_max', 'moments_hu', 'perimeter',
           'perimeter_crofton', 'solidity', 'bbox', 'centroid')


def curvature(this_mask, Lval):
    try:
        outline = (skimage.morphology.dilation(this_mask) ^ this_mask)*1.0
        curves = skimage.measure.find_contours(outline, 0.0, fully_connected='high')
        # Smooth contour with Gaussian kernel
        sigma = 3.0

        smoothedCurves={}
        for c in range(0, len(curves)):
            sx = nd.gaussian_filter1d(curves[c][:,0], sigma)
            sy = nd.gaussian_filter1d(curves[c][:,1], sigma)

            smoothedCurves[c] = np.vstack((sx,sy))
        objectSizes=[len(smoothedCurves[c][0]) for c in smoothedCurves]
        ind=np.where(objectSizes==np.max(objectSizes))[0][0]

        # Input parameters for CPDA
    #    Lvals=[10,25,50,75,100,200,400,800,1000]
        Lvals=[Lval]
        # Compute CPDA
        for L in Lvals:
            # If indexing runs clockwise, flip it:
            curve=np.flip(smoothedCurves[ind],axis=1)
            # Else, no need to do anything:
            # curve=smoothedCurves[ind]
            H_L=cpda.cpda(curve, L) # This computes curvature

        mean = np.mean(H_L)
        std = np.std(H_L)
        mymin = np.min(H_L)
        mymax = np.max(H_L)
        frac_below_zero = np.sum(H_L < 0)/len(H_L)
        

        return mean, std, mymin, mymax, frac_below_zero
    except:
        return -1.0, -1.0, -1.0, -1.0, -1.0
    
def calculate_cpmeasure(pixels, masks):
    """
    Calculate cell morphology measurements using cp_measure.

    Parameters:
    - pixels: 2D numpy array representing the image pixels.
    - masks: 2D numpy array representing the labeled masks.

    Returns:
    - results: DataFrame containing the morphology measurements for each labeled region.
    """

    masks, fw, inv = skimage.segmentation.relabel_sequential(masks)
    measurements = get_core_measurements()

    results = {}
    for name, v in measurements.items():
        results = {**results, **v(masks, pixels)}

    df_results = pd.DataFrame(results)
    df_results.index = df_results.index + 1
    df_results.index = df_results.index.map(inv)
    return df_results

def is_on_edge2(this_image, df2, this_cellID):
    frame_shape0 = df2.loc[this_cellID, 'frame_shape0'].values[0]
    frame_shape1 = df2.loc[this_cellID, 'frame_shape1'].values[0]
    bbox0 = df2.loc[this_cellID, 'BoundingBoxMinimum_X'].values[0]
    bbox1 = df2.loc[this_cellID, 'BoundingBoxMinimum_Y'].values[0]
    bbox2 = df2.loc[this_cellID, 'BoundingBoxMaximum_X'].values[0]
    bbox3 = df2.loc[this_cellID, 'BoundingBoxMaximum_Y'].values[0]
    
    margin = 20
    if bbox0 < margin:
        test_area = this_image[0:bbox0, bbox1:bbox3]
        if np.all(test_area == 0):
            return 1
    if bbox1 < margin:
        test_area = this_image[bbox0:bbox2, 0:bbox1]
        if np.all(test_area == 0):
            return 1
    if frame_shape0 - bbox2 < margin:
        test_area = this_image[bbox2:, bbox1:bbox3]
        if np.all(test_area == 0):
            return 1
    if frame_shape1 - bbox3 < margin:
        test_area = this_image[bbox0:bbox2, bbox3:]
        if np.all(test_area == 0):
            return 1
    return 0

list_of_dataframes = []
Lvals = [5, 10, 20, 50, 100, 200]
curvature_funcs = [partial(curvature, Lval=Lval) for Lval in Lvals]

# Add __name__ attribute to partial objects for skimage compatibility
for Lval, func in zip(Lvals, curvature_funcs):
    func.__name__ = f'curvature-{Lval}'

# This part of the code iterates over frames to capture cell information
for this_frame_id in unique_frame_ids:
    try:
        this_frame = seg_relabeled[int(this_frame_id)]
        this_image = images[int(this_frame_id)]
    except:
        print("Error in frame ", this_frame_id)
        continue
    
    df = pd.DataFrame(regionprops_table(this_frame, properties = ('label',), extra_properties=curvature_funcs))
    df.index = df['label']
    df2 = calculate_cpmeasure(this_image, this_frame)
    df = pd.merge(df, df2, left_index=True, right_index=True)
    list_of_cellIDs = list(df['label'])
    df['frame_id_T'] = int(this_frame_id)
    df['realTime'] = int(this_frame_id) * cycleTime
    df['frame_shape0'] = this_frame.shape[0]
    df['frame_shape1'] = this_frame.shape[1]
    df['Difference_of_Location_CenterMassIntensity_X_to_BoundingBoxMinimum_X'] =\
        df['Location_CenterMassIntensity_X'] - df['BoundingBoxMinimum_X']
    df['Difference_of_Location_CenterMassIntensity_Y_to_BoundingBoxMinimum_Y'] =\
        df['Location_CenterMassIntensity_Y'] - df['BoundingBoxMinimum_Y']
    df['Difference_of_Location_MaxIntensity_X_to_BoundingBoxMinimum_X'] =\
        df['Location_MaxIntensity_X'] - df['BoundingBoxMinimum_X']
    df['Difference_of_Location_MaxIntensity_Y_to_BoundingBoxMinimum_Y'] =\
        df['Location_MaxIntensity_Y'] - df['BoundingBoxMinimum_Y']
    df['Difference_of_Center_X_to_BoundingBoxMinimum_X'] =\
        df['Center_X'] - df['BoundingBoxMinimum_X']
    df['Difference_of_Center_Y_to_BoundingBoxMinimum_Y'] =\
        df['Center_Y'] - df['BoundingBoxMinimum_Y']
    df['Length_of_BoundingBox_X'] =\
        df['BoundingBoxMaximum_X'] - df['BoundingBoxMinimum_X']
    df['Length_of_BoundingBox_Y'] =\
        df['BoundingBoxMaximum_Y'] - df['BoundingBoxMinimum_Y']


    rename_dict = {'label':'cellID'}
    for j, stat in enumerate(['mean', 'std', 'min', 'max', 'frac_below_zero']):
        rename_dict.update({f'curvature-{Lval}-{j}': f'curvature_{Lval}-{stat}' for Lval in Lvals})

    df2 = df.rename(rename_dict, axis=1).set_index(['cellID', 'frame_id_T'])
    
#    mycellid_props = df2.to_dict('index')
#    mygraph = skimage.graph.RAG(this_frame)
#    mygraph.remove_node(0)

    this_frame_data = data[data[:,1] == this_frame_id]
    for this_cellID in list_of_cellIDs:
        this_cellID_frame_data = this_frame_data[this_frame_data[:,0] == this_cellID]
        if len(this_cellID_frame_data) == 1:
            df2.loc[this_cellID, 'track_pos0'] = this_cellID_frame_data[0, 2]
            df2.loc[this_cellID, 'track_pos1'] = this_cellID_frame_data[0, 3]
        else:
            df2.loc[this_cellID, 'track_pos0'] = -1
            df2.loc[this_cellID, 'track_pos1'] = -1
        
        df2.loc[this_cellID, 'number_1px_away'] = int(get_number_of_neighbors(this_cellID, this_frame, pixel_range=1))
        df2.loc[this_cellID, 'is_on_edge2'] = is_on_edge2(this_image, df2, this_cellID)
#        df2.loc[this_cellID, 'number_of_adjacent_objects'] = mygraph.degree(this_cellID)

    list_of_dataframes.append(df2)

df_all = pd.concat(list_of_dataframes)

def get_instantaneous_speeds(this_cell_df, cycleTime, step):
    """
    This function calculates the instantaneous speed of the cell
    The speed is calculated as the distance between the next position and the current position
    divided by the real time difference between the next frame and the current frame
    This_cell_df is a pandas dataframe that contains the information of one cell through time
    It must have columns for Center_X and Center_Y
    The index should be the index of time to be multiplied by cycleTime to get the real time
    This function does not currently handle the resolution of the image,
    but returns the speed in pixels per real time.

    This functions leaves the last values as -1 because there is no next positions to calculate the speed
    """
    
    x_positions = this_cell_df['Center_X'].values
    y_positions = this_cell_df['Center_Y'].values
    times = (this_cell_df.index*cycleTime).values
    speed = -1*np.ones(len(x_positions), float)
    speed_angle = -1*np.ones(len(x_positions), float)
    try:
        for i in range(len(x_positions) -1):
            displacement = (x_positions[i+step] - x_positions[i], 
                            y_positions[i+step] - y_positions[i])
            dist = np.linalg.norm(displacement)
            time_difference = times[i+step] - times[i]
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


def measure_directionality_ratio(this_cell_df, cycleTime):
    """
    This function calculates the directionality ratio of the cell. 
    This is calculated as the ratio of the distance from the origin to the current position 
    divided by the total distance traversed by the cell (the path distance).


    """


    # This function sets the first value as 1 and the second value is 1 by calculation
    x_positions = this_cell_df['Center_X'].values
    y_positions = this_cell_df['Center_Y'].values
    start_position = (x_positions[0], y_positions[0])
    total_distance_traversed = 0
    directionality_ratios = -1 * np.ones(len(x_positions), float)
    try:
        for i in range(len(x_positions)):
            if i == 0:
                directionality_ratios[i] = 1
            else:
                displacement_since_last_position = (x_positions[i] - x_positions[i-1], 
                                                    y_positions[i] - y_positions[i-1])
                dist_since_last_position = np.linalg.norm(displacement_since_last_position)
                total_distance_traversed += dist_since_last_position

                displacement_since_origin = (x_positions[i] - x_positions[0], 
                                                    y_positions[i] - y_positions[0])
                dist_since_origin = np.linalg.norm(displacement_since_origin)

                directionality_ratios[i] = dist_since_origin / total_distance_traversed
        return directionality_ratios
    except:
        return directionality_ratios


def get_displacement_norm_from_origin(this_cell_df):
    """"
    This function calculates the displacement of the cell from the origin
    It ignores the intermediate path
    """

    x_positions = this_cell_df['Center_X'].values
    y_positions = this_cell_df['Center_Y'].values
    start_position = (x_positions[0], y_positions[0])
    displacement_norm_from_origin = -1 * np.ones(len(x_positions), float)
    displacement_norm_from_origin_squared = -1 * np.ones(len(x_positions), float)
    try:
        for i in range(len(x_positions)):
            if i == 0:
                displacement_norm_from_origin[i] = 0
                displacement_norm_from_origin_squared[i] = 0
            else:
                displacement_since_origin = (x_positions[i] - x_positions[0], 
                                                            y_positions[i] - y_positions[0])
                this_displacement_norm_since_origin = np.linalg.norm(displacement_since_origin)
                displacement_norm_from_origin[i] = this_displacement_norm_since_origin
                displacement_norm_from_origin_squared[i] = this_displacement_norm_since_origin**2
        return displacement_norm_from_origin, displacement_norm_from_origin_squared
    except:
        return displacement_norm_from_origin, displacement_norm_from_origin_squared


def get_path_distance_from_origin(this_cell_df):
    """
    This function calculates the path distance of the cell from the origin
    This is the total distance traversed by the cell
    """

    x_positions = this_cell_df['Center_X'].values
    y_positions = this_cell_df['Center_Y'].values
    start_position = (x_positions[0], y_positions[0])
    total_distance_traversed = 0
    total_distance_traversed_array = -1 * np.ones(len(x_positions), float)
    try:
        for i in range(len(x_positions)):
            if i == 0:
                total_distance_traversed_array[i] = 0
            else:
                displacement_since_last_position = (x_positions[i] - x_positions[i-1], 
                                                    y_positions[i] - y_positions[i-1])
                dist_since_last_position = np.linalg.norm(displacement_since_last_position)
                total_distance_traversed += dist_since_last_position
                total_distance_traversed_array[i] = total_distance_traversed
                
        return total_distance_traversed_array
    except:
        return total_distance_traversed_array
    

cell_indices = list(df_all.index.get_level_values(0).unique())

#This part of the code iterates over cells to calculate the dynamics features
for this_cellID in cell_indices:
#    if this_cellID != 3:
#        continue
    this_cell_df = df_all.loc[this_cellID]
    for i in range(5):
        step = i + 1
        speed, speed_angle = get_instantaneous_speeds(this_cell_df, cycleTime, step)
        df_all.loc[this_cellID, f'speed_step_{step}'] = speed
        df_all.loc[this_cellID, f'speed_angle_step_{step}'] = speed_angle
    df_all.loc[this_cellID, 'N_frames_existence'] = int(len(this_cell_df))
    df_all.loc[this_cellID, 'directionality_ratio'] = measure_directionality_ratio(this_cell_df, cycleTime)
    displacement_norm_from_origin, displacement_norm_from_origin_squared =\
    get_displacement_norm_from_origin(this_cell_df)
    df_all.loc[this_cellID, 'displacement_norm_from_origin'] = displacement_norm_from_origin
    df_all.loc[this_cellID, 'displacement_norm_from_origin_squared'] = displacement_norm_from_origin_squared
    df_all.loc[this_cellID, 'path_distance_from_origin'] = get_path_distance_from_origin(this_cell_df)
#    display(this_cell_df)
#display(df_all)

first_cols = ['is_on_edge2', 'frame_shape0', 'frame_shape1',
                'Location_CenterMassIntensity_X', 'Location_CenterMassIntensity_Y', 
                'Location_CenterMassIntensity_Z', 'Location_MaxIntensity_X', 
                'Location_MaxIntensity_Y', 'Location_MaxIntensity_Z',
                'BoundingBoxMinimum_X', 'BoundingBoxMaximum_X', 'BoundingBoxMinimum_Y', 'BoundingBoxMaximum_Y',
                'track_pos0', 'track_pos1',   'Center_X', 'Center_Y', 'N_frames_existence', 'number_1px_away', 'realTime',
                'EulerNumber',]
other_cols = [col for col in df_all.columns if col not in first_cols]
df_all = df_all[first_cols + other_cols]
df_all.to_csv(this_output)