import os
import natsort
import re
files = glob_wildcards("1_data/{subfolder_filename}.tif")
both = glob_wildcards("1_data/{subfolder}/{filename}.tif")
both.subfolder
both.filename

rule all:
    input:
#        expand("3c_tracking_images_filtered/{subfolder_filename}.tif", subfolder_filename = files.subfolder_filename)
        expand("4b_time_averaged_cell_morphodynamics/{subfolder}/cell_data.csv", subfolder = both.subfolder),
        expand("5_tracking_images_outlines/{subfolder_filename}.tif", subfolder_filename = files.subfolder_filename)

####################################################################################################
# Preprocessing

def get_stabilization_files_list_from_subfolder(wildcards):
    this_original_sub = "1_data/" + wildcards.subfolder
    this_root_sub = "1_data/" + wildcards.subfolder
    list_of_stabilized_images = [os.path.join(this_root_sub, each) for 
        each in natsort.natsorted(os.listdir(this_original_sub))]
    return list_of_stabilized_images

#rule stabilize_individuals_to_group:
#    input:
#        get_stabilization_files_list_from_subfolder
#    output:
#        "1b_stabilize/{subfolder}/{filename}.tif"
#        ["1b_stabilize/{subfolder}/{filename}.txt".format(filename=filename) for filename in os.listdir("1b_stabilize/{subfolder}/")]
#        "1b_stabilize/{subfolder}/image_stack.npy"
#    script:
#        "scripts/run_pystackreg.py"

def get_image_stack_from_subfolderfilename(wildcards):
    subfolder = os.path.split(os.path.split(wildcards.subfolder_filename)[0])[1]
    image_stack = os.path.join('1b_stabilize', subfolder, 'image_stack.npy')
    return image_stack

#rule stabilize_group_back_to_individuals:
#    input:
#        get_image_stack_from_subfolderfilename,
#        "1_data/{subfolder_filename}.tif"
#    output:
#        "1c_stabilize/{subfolder_filename}.tif"
#    script:
#        "scripts/run_group_back_to_individuals.py"
    

####################################################################################################
# Segmentation

#rule segment_with_stardist:
#    input:
#        "1_data/{subfolder_filename}.tif"
#    output:
#        "2_segmentation/{subfolder_filename}.tif"
#    retries: 10
#    shell:
#        "python scripts/run_stardist.py {input} {output}"


def get_equivalent_nuclear_segmentation(wildcards):
#    this_root_sub = os.path.abspath("../../2024-03-14_snakemake_develop_tracking/MorphoDynamicsPipe/2_segmentation/")
    this_root_sub = os.path.abspath("../../2024-04-25_new_lifs_batch_nucs/MorphoDynamicsPipe/1_data/")
    new_subfolder_filename = re.sub('C=2', 'C=0', wildcards.subfolder_filename)
    new_fullpath = os.path.join(this_root_sub, new_subfolder_filename) + '.tif'
#    print(new_fullpath)
    return [new_fullpath]

#rule segment_with_micro_sam:
#    input:
#        "1_data/{subfolder_filename}.tif",
#        get_equivalent_nuclear_segmentation
#    output:
#        "2_segmentation/{subfolder_filename}.tif"
#    retries: 10
#    script:
#        "scripts/run_microsam.py"

rule segment_with_cellpose_nucs:
    input:
        "1_data/{subfolder_filename}.tif",
    output:
        "2_segmentation/{subfolder_filename}.tif"
    script:
        "scripts/run_cellpose_nucs.py"

#rule segment_with_cellpose_celltracker_with_nucs:
#    input:
#        "1_data/{subfolder_filename}.tif",
#        get_equivalent_nuclear_segmentation #this is the link to the nuclear channel
#    output:
#        "2_segmentation/{subfolder_filename}.tif"
#    script:
#        "scripts/run_cellpose_celltracker_with_nucs.py"

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

def get_tracking_info_from_subfolderfilename_nuclei_version(wildcards):
    subfolder = os.path.split(os.path.split(wildcards.subfolder_filename)[0])[1]
    subfolder_renamed = re.sub('C=2', 'C=0', subfolder)
    tracking_info_base_path = os.path.abspath('../../2024-04-25_new_lifs_batch_nucs/MorphoDynamicsPipe/3a_tracking_info')
#    tracking_info_base_path = os.path.abspath('3a_tracking_info')
    tracking_info = os.path.join(tracking_info_base_path, subfolder_renamed, 'track_info.npy')
    return tracking_info

rule convert_btrack_info_to_images:
    input:
        "2_segmentation/{subfolder_filename}.tif",
        get_tracking_info_from_subfolderfilename,
    output:
#        "3b_tracking_images/{subfolder_filename}.tif"
        "3b_tracking_images_from_nucs/{subfolder_filename}.tif"
    script:
        "scripts/convert_btrack_info_to_images.py"

rule filter_cells_after_tracking:
    input:
#        "3b_tracking_images/{subfolder_filename}.tif",
        "3b_tracking_images_from_nucs/{subfolder_filename}.tif",
        get_tracking_info_from_subfolderfilename,
    output:
        "3c_tracking_images_filtered/{subfolder_filename}.tif"
    script:
        "scripts/filter_short_lived_cells_after_tracking.py"


####################################################################################################
# Cell Morphodynamics

def get_segmentation_relabeled_files_list_from_subfolder(wildcards):
    this_original_sub = "1_data/" + wildcards.subfolder
#    this_root_sub = "3b_tracking_images/" + wildcards.subfolder
    this_root_sub = "3b_tracking_images_from_nucs/" + wildcards.subfolder
    list_of_segmented_images = [os.path.join(this_root_sub, each)
        for each in natsort.natsorted(os.listdir(this_original_sub))]
    return list_of_segmented_images

rule extract_instantaneous_cell_morphodynamics:
    input:
        "3a_tracking_info/{subfolder}/track_info.npy",
        get_segmentation_relabeled_files_list_from_subfolder,
    output:
        "4a_instantaneous_cell_morphodynamics/{subfolder}/cell_data.csv"
    retries: 2
    script:
        "scripts/extract_instantaneous_cell_morphodynamics.py"

def get_list_of_input_subfolders(wildcards):
    list_of_subfolders = [os.path.join('1_data', each) for each in os.listdir('1_data') if os.path.isdir(os.path.join('1_data', each))]
    return list_of_subfolders

rule get_maximum_common_time:
    input:
        get_list_of_input_subfolders
    output:
        "maximum_common_time.txt"
    script:
        "scripts/get_maximum_common_time.py"


rule extract_time_averaged_cell_morphodynamics:
    input:
        "4a_instantaneous_cell_morphodynamics/{subfolder}/cell_data.csv",
        "maximum_common_time.txt",
        get_segmentation_relabeled_files_list_from_subfolder,
    output:
        "4b_time_averaged_cell_morphodynamics/{subfolder}/cell_data.csv"
    script:
        "scripts/extract_time_averaged_cell_morphodynamics.py"

####################################################################################################
# Visualization (optional)

rule convert_to_label_outlines:
    input:
        "3b_tracking_images_from_nucs/{subfolder_filename}.tif"
    output:
        '5_tracking_images_outlines/{subfolder_filename}.tif'
    script:
        "scripts/convert_to_label_outlines.py"
