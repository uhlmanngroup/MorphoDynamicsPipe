import os
import natsort
files = glob_wildcards("1_data/{subfolder_filename}.tif")
both = glob_wildcards("1_data/{subfolder}/{filename}.tif")
both.subfolder
both.filename

rule all:
    input:
#        expand("3b_tracking_images/{subfolder_filename}.tif", subfolder_filename = files.subfolder_filename)
        expand("4b_time_averaged_cell_morphodynamics/{subfolder}/cell_data.csv", subfolder = both.subfolder)

####################################################################################################
# Segmentation

rule segment_with_stardist:
    input:
        "1_data/{subfolder_filename}.tif"
    output:
        "2_segmentation/{subfolder_filename}.tif"
    retries: 10
    shell:
        "python scripts/run_stardist.py {input} {output}"

####################################################################################################
# Tracking

def get_segmentation_files_list_from_subfolder(wildcards):
    this_original_sub = "1_data/" + wildcards.subfolder
    this_root_sub = "2_segmentation/" + wildcards.subfolder
    list_of_segmented_images = [os.path.join(this_root_sub, each) for 
        each in natsort.natsorted(os.listdir(this_original_sub))]
    return list_of_segmented_images


rule track_with_btrack:
    input:
        get_segmentation_files_list_from_subfolder
    output:
        "3a_tracking_info/{subfolder}/track_info.npy"
    script:
        "scripts/run_btrack_to_info.py"


def get_tracking_info_from_subfolderfilename(wildcards):
    subfolder = os.path.split(os.path.split(wildcards.subfolder_filename)[0])[1]
    tracking_info = os.path.join('3a_tracking_info', subfolder, 'track_info.npy')
    return tracking_info

rule convert_btrack_info_to_images:
    input:
        "2_segmentation/{subfolder_filename}.tif",
        get_tracking_info_from_subfolderfilename,
    output:
        "3b_tracking_images/{subfolder_filename}.tif"
    script:
        "scripts/convert_btrack_info_to_images.py"

####################################################################################################
# Cell Morphodynamics

def get_segmentation_relabeled_files_list_from_subfolder(wildcards):
    this_original_sub = "1_data/" + wildcards.subfolder
    this_root_sub = "3b_tracking_images/" + wildcards.subfolder
    list_of_segmented_images = [os.path.join(this_root_sub, each)
        for each in natsort.natsorted(os.listdir(this_original_sub))]
    return list_of_segmented_images

rule extract_instantaneous_cell_morphodynamics:
    input:
        "3a_tracking_info/{subfolder}/track_info.npy",
        get_segmentation_relabeled_files_list_from_subfolder,
    output:
        "4a_instantaneous_cell_morphodynamics/{subfolder}/cell_data.csv"
    script:
        "scripts/extract_instantaneous_cell_morphodynamics.py"


rule extract_time_averaged_cell_morphodynamics:
    input:
        "4a_instantaneous_cell_morphodynamics/{subfolder}/cell_data.csv",
        get_segmentation_relabeled_files_list_from_subfolder,
    output:
        "4b_time_averaged_cell_morphodynamics/{subfolder}/cell_data.csv"
    script:
        "scripts/extract_time_averaged_cell_morphodynamics.py"