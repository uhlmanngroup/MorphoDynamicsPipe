# This is a snakemake file to run the MorphoDynamicsPipe pipeline using the example data provided in the repository.
# The example images must be copied to 1_data folder before running this pipeline.
# uses morphody50 conda environment
# snakemake -s run_example.smk --cores "all" --keep-going
# snakemake -s run_example.smk --cores 4 --keep-going
# or on slurm
# sbatch -t 24:00:00 --mem=64G -c 16 --gres=gpu:1 --wrap="snakemake -s run_example.smk --cores "all" --keep-going"

import os
import natsort
import re
import platform
import subprocess
print('My platform is: ', platform.system())

####################################################################################################
# Defining files to create
files_wcs = glob_wildcards("1_data/{subfolder_filename}.tif")
# filtered list (do not overwrite files_wcs attributes)
files_subfolder_filenames = [f for f in files_wcs.subfolder_filename if not f.endswith(":Zone.Identifier")]

both_wcs = glob_wildcards("1_data/{subfolder}/{filename}.tif")
both_subfolders = [f for f in both_wcs.subfolder if not f.endswith(":Zone.Identifier")]
both_filenames = [f for f in both_wcs.filename if not f.endswith(":Zone.Identifier")]
# (use files_subfolder_filenames, both_subfolders, both_filenames below)

rule all:
    input:
    #commenting these lines in and out will control how many steps the pipeline performs
        expand("2_segmentation/{subfolder_filename}.tif", subfolder_filename = files_subfolder_filenames),
        expand("3b_tracking_images/{subfolder_filename}.tif", subfolder_filename = files_subfolder_filenames),
        expand("3c_tracking_images_filtered/{subfolder_filename}.tif", subfolder_filename = files_subfolder_filenames),
        expand("4a_instantaneous_cell_morphodynamics/{subfolder}/cell_data.csv", subfolder = both_subfolders),
        expand("4b_time_averaged_cell_morphodynamics/{subfolder}/cell_data.csv", subfolder = both_subfolders),
        expand("5_tracking_images_outlines/{subfolder_filename}.tif", subfolder_filename = files_subfolder_filenames)

####################################################################################################
# Preprocessing

def get_stabilization_files_list_from_subfolder(wildcards):
    this_original_sub = "1_data/" + wildcards.subfolder
    this_root_sub = "1_data/" + wildcards.subfolder
    list_of_stabilized_images = [this_root_sub + '/' + each for 
        each in natsort.natsorted([each2 for each2 in os.listdir(this_original_sub) if not each2.startswith('.')])]
    return list_of_stabilized_images

def get_image_stack_from_subfolderfilename(wildcards):
    subfolder = os.path.split(os.path.split(wildcards.subfolder_filename)[0])[1]
    image_stack = '1b_stabilize/', subfolder, '/image_stack.npy'
    return image_stack

#####

#rule stabilize_individuals_to_group:
#    input:
#        get_stabilization_files_list_from_subfolder
#    output:
#        "1b_stabilize/{subfolder}/{filename}.tif"
#        ["1b_stabilize/{subfolder}/{filename}.txt".format(filename=filename) for filename in [each2 for each2 in os.listdir("1b_stabilize/{subfolder}/") if not each2.startswith('.')]]
#        "1b_stabilize/{subfolder}/image_stack.npy"
#    script:
#        "scripts/run_pystackreg.py"

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

def get_equivalent_nuclear_segmentation(wildcards):
#    this_root_sub = os.path.abspath("../../2024-03-14_snakemake_develop_tracking/MorphoDynamicsPipe/2_segmentation/")
    this_root_sub = os.path.abspath("../../2024-07-30_new_lifs_batch_nucs_slightly_broader_dataset/MorphoDynamicsPipe/1_data/")
    new_subfolder_filename = re.sub('C=2', 'C=0', wildcards.subfolder_filename)
    new_fullpath = os.path.join(this_root_sub, new_subfolder_filename) + '.tif'
    new_fullpath = new_fullpath.replace(os.sep, '/')
#    print(new_fullpath)
    return [new_fullpath]

#####
# comment in exactly one of the following segmentation rules

# comment this in to run cellpose on a single channel, eg. celltracker alone
rule segment_with_cellpose_nucs:
    input:
        "1_data/{subfolder_filename}.tif",
    output:
        "2_segmentation/{subfolder_filename}.tif"
    retries: 0
    shell:
        "python scripts/run_cellpose_nucs_v4.py {input} {output}"


####################################################################################################
# Tracking

def get_segmentation_files_list_from_subfolder(wildcards):
    this_original_sub = "1_data/" + wildcards.subfolder
    this_root_sub = "2_segmentation/" + wildcards.subfolder
    list_of_segmented_images = [this_root_sub + '/' + each for
        each in natsort.natsorted([each2 for each2 in os.listdir(this_original_sub) if not each2.startswith('.')])]
    return list_of_segmented_images

def get_tracking_info_from_subfolderfilename(wildcards):
    subfolder = os.path.split(os.path.split(wildcards.subfolder_filename)[0])[1]
    tracking_info = '3a_tracking_info/' + subfolder + '/track_info.npy'
    return tracking_info

def get_tracking_info_from_subfolderfilename_nuclei_version(wildcards):
    subfolder = os.path.split(os.path.split(wildcards.subfolder_filename)[0])[1]
    subfolder_renamed = re.sub('C=2', 'C=0', subfolder)
    tracking_info_base_path = os.path.abspath('../../2024-07-30_new_lifs_batch_nucs_slightly_broader_dataset/MorphoDynamicsPipe/3a_tracking_info')
#    tracking_info_base_path = os.path.abspath('3a_tracking_info')
    tracking_info = os.path.join(tracking_info_base_path, subfolder_renamed, 'track_info.npy')
    tracking_info = tracking_info.replace(os.sep, '/')
    return tracking_info

def get_tracking_info_from_subfolder_nuclei_version(wildcards):
    subfolder = wildcards.subfolder
    subfolder_renamed = re.sub('C=2', 'C=0', subfolder)
    tracking_info_base_path = os.path.abspath('../../2024-07-30_new_lifs_batch_nucs_slightly_broader_dataset/MorphoDynamicsPipe/3a_tracking_info')
#    tracking_info_base_path = os.path.abspath('3a_tracking_info')
    tracking_info = os.path.join(tracking_info_base_path, subfolder_renamed, 'track_info.npy')
    tracking_info = tracking_info.replace(os.sep, '/')
    return tracking_info

#####

# comment in for nuclei or celltracker without nuclei and out for celltracker with nuclei
rule track_with_btrack:
    input:
        get_segmentation_files_list_from_subfolder
    output:
        "3a_tracking_info/{subfolder}/track_info.npy"
#    container:
#        "docker://spectralnanodiamond/btrack:latest"
#        "../../2024-08-22_making_example/MorphoDynamicsPipe/.snakemake/singularity/443cbde37592944ef7c547806b1792f4.simg"
#    conda:
#        btrack_conda_env
    script:
        "scripts/run_btrack_to_info.py"

# comment in for celltracker version that relies on nuclei and out otherwise
#rule symlink_to_btrack_info:
#choose one of the shell scripts to run depending on your operating system
#    input:
#        get_tracking_info_from_subfolder_nuclei_version
#    output:
#        "3a_tracking_info/{subfolder}/track_info.npy"
#    shell:
#        "ln -rs {input} {output}" #mac / linux
#        "cp {input} {output}"    #alternative mac / linux
#        "python -c 'import os; os.symlink({input}, {output})'" #this hopefully works on windows

rule convert_btrack_info_to_images:
    input:
        "2_segmentation/{subfolder_filename}.tif",
        get_tracking_info_from_subfolderfilename
    output:
        "3b_tracking_images/{subfolder_filename}.tif"
    script:
        "scripts/convert_btrack_info_to_images.py"

rule filter_cells_after_tracking:
    input:
        "3b_tracking_images/{subfolder_filename}.tif",
        get_tracking_info_from_subfolderfilename,
    output:
        "3c_tracking_images_filtered/{subfolder_filename}.tif"
    script:
        "scripts/filter_short_lived_cells_after_tracking.py"


#use this rule to create trivial tracking info without time course data
#rule create_trivial_btrack_info:
#    input:
#        get_segmentation_files_list_from_subfolder
#    output:
#        "3a_tracking_info/{subfolder}/track_info.npy"
#    container:
#        "docker://spectralnanodiamond/btrack:latest"
#        "../../2024-08-22_making_example/MorphoDynamicsPipe/.snakemake/singularity/443cbde37592944ef7c547806b1792f4.simg"
#    script:
#        "scripts/create_trivial_btrack_info.py"

####################################################################################################
# Cell Morphodynamics

def get_segmentation_relabeled_files_list_from_subfolder(wildcards):
    this_original_sub = "1_data/" + wildcards.subfolder
#    this_root_sub = "3b_tracking_images/" + wildcards.subfolder
    this_root_sub = "3c_tracking_images_filtered/" + wildcards.subfolder
    list_of_segmented_images = [this_root_sub + '/' + each
        for each in natsort.natsorted([each2 for each2 in os.listdir(this_original_sub) if not each2.startswith('.')])]
    return list_of_segmented_images

def get_images_files_list_from_subfolder(wildcards):
    this_original_sub = "1_data/" + wildcards.subfolder
    list_of_images = [this_original_sub + '/' + each
        for each in natsort.natsorted([each2 for each2 in os.listdir(this_original_sub) if not each2.startswith('.')])]
    return list_of_images

def get_list_of_input_subfolders(wildcards):
    list_of_subfolders = ['1_data/' + each for each in os.listdir('1_data') if os.path.isdir(os.path.join('1_data', each))]
    return list_of_subfolders

#####

rule extract_instantaneous_cell_morphodynamics:
    input:
        "3a_tracking_info/{subfolder}/track_info.npy",
        get_segmentation_relabeled_files_list_from_subfolder,
        get_images_files_list_from_subfolder,
    output:
        "4a_instantaneous_cell_morphodynamics/{subfolder}/cell_data.csv"
    retries: 2
    script:
        "scripts/extract_instantaneous_cell_morphodynamics.py"

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
#        "3b_tracking_images/{subfolder_filename}.tif"
        "3c_tracking_images_filtered/{subfolder_filename}.tif"
    output:
        '5_tracking_images_outlines/{subfolder_filename}.tif'
    script:
        "scripts/convert_to_label_outlines.py"
