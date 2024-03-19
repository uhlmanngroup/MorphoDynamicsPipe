import os
import natsort
files = glob_wildcards("1_data/{subfolder_filename}.tif")
both = glob_wildcards("1_data/{subfolder}/{filename}.tif")
both.subfolder
both.filename

rule all:
    input:
#        expand("3b_tracking_images/{subfolder_filename}.tif", subfolder_filename = files.subfolder_filename)
        expand("4_cell_morphodynamics/{subfolder}/cell_data.csv", subfolder = both.subfolder)

rule segment_with_stardist:
    input:
        "1_data/{subfolder_filename}.tif"
    output:
        "2_segmentation/{subfolder_filename}.tif"
    retries: 10
    shell:
        "python scripts/run_stardist.py {input} {output}"


def get_subfolder_files_list(wildcards):
    this_original_sub = "1_data/" + wildcards.subfolder
    this_root_sub = "2_segmentation/" + wildcards.subfolder
    list_of_segmented_images = [this_root_sub + '/' + each for each in natsort.natsorted(os.listdir(this_original_sub))]
    print('this is in the snakemake function 1')
    print(list_of_segmented_images)
    return list_of_segmented_images

rule track_with_btrack:
    input:
        get_subfolder_files_list
    output:
        "3a_tracking_info/{subfolder}/track_info.npy"
    script:
        "scripts/run_btrack_to_info.py"


def get_tracking_info_from_subfolderfilename(wildcards):
    subfolder = os.path.split(os.path.split(wildcards.subfolder_filename)[0])[1]
    tracking_info = '3a_tracking_info/' + subfolder + '/track_info.npy'
    print('this is in the snakemake function 2')
    print(tracking_info)
    return tracking_info

rule convert_btrack_info_to_images:
    input:
        "2_segmentation/{subfolder_filename}.tif",
        get_tracking_info_from_subfolderfilename,
    output:
        "3b_tracking_images/{subfolder_filename}.tif"
    script:
        "scripts/convert_btrack_info_to_images.py"


def get_segmentation_relabeled_list(wildcards):
    this_original_sub = "1_data/" + wildcards.subfolder
    this_root_sub = "3b_tracking_images/" + wildcards.subfolder
    list_of_segmented_images = [this_root_sub + '/' + each for each in natsort.natsorted(os.listdir(this_original_sub))]
    print('this is in the snakemake function 3')
    print(list_of_segmented_images)
    return list_of_segmented_images

rule extract_cell_morphodynamics:
    input:
        "3a_tracking_info/{subfolder}/track_info.npy",
        get_segmentation_relabeled_list,
    output:
        "4_cell_morphodynamics/{subfolder}/cell_data.csv"
    script:
        "scripts/extract_cell_morphodynamics.py"


#rule get_regionprops:
#    input:
#        "2_segmentation/{subfolder_filename}.tif"
#    output:
#        "4_regionprops/{subfolder_filename}.csv"
#    conda:
#        "conda_envs_yaml/environment_StardistSmake_dev.yml"
#    shell:
#        "python scripts/get_regionprops.py {input} {output}"

#    shell:
#        "python scripts/multiple_inputs.py {input} {output}"
# I think you could also use shell here if you parse the list thing into a space separated string or something