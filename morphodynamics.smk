files = glob_wildcards("1_data/{subfolder_filename}.tif")
both = glob_wildcards("1_data/{subfolder}/{filename}.tif")
both.subfolder
both.filename

rule all:
    input:
        expand("3b_tracking_images/{subfolder_filename}.tif", subfolder_filename = files.subfolder_filename)

#rule segment_with_stardist:
#    input:
#        "1_data/{subfolder_filename}.tif"
#    output:
#        "2_segmentation/{subfolder_filename}.tif"
#    retries: 10
#    conda:
#        "conda_envs_yaml/environment_StardistSmake_dev.yml"
#    shell:
#        "python scripts/run_stardist.py {input} {output}"

#rule track_with_btrack:
#    input:
#        "../data/{subfolder}
#        expand("../data/{subfolder}/{{file}}.tif", subfolder = both.subfolder)
# double curly brackets use by different rules
#    output:
#        expand("results/{subfolder_filename}.tif", subfold
# if output is directory only then it needs to be wrapped in directory only
#    conda:
#        "conda_envs_yaml/environment_BtrackSmake_dev.yml"
#    script:
#        scripts/run_btrack.py
# snakemake.input
#

#rule get_regionprops:
#    input:
#        "2_segmentation/{subfolder_filename}.tif"
#    output:
#        "4_regionprops/{subfolder_filename}.csv"
#    conda:
#        "conda_envs_yaml/environment_StardistSmake_dev.yml"
#    shell:
#        "python scripts/get_regionprops.py {input} {output}"

def get_subfolder_files_list(wildcards):
#    this = os.listdir("1_data/" + both.subfolder)
    this_root_sub = "1_data/" + wildcards.subfolder
    this2 = [this_root_sub + '/' + each for each in os.listdir(this_root_sub)]
    print('this is in the snakemake function 1')
    print(this2)
    return this2

rule test_multiple_inputs_to_one_output:
    input:
        get_subfolder_files_list
#        os.listdir("1_data/{subfolder}".format(subfolder = subfolder) for subfolder in both.subfolder)
    output:
        "3a_tracking_info/{subfolder}/output.npy"
#    conda:
#        "conda_envs_yaml/environment_StardistSmake_dev.yml"
#    shell:
#        "python scripts/multiple_inputs.py {input} {output}"
# I think you could also use shell here if you parse the list thing into a space separated string or something
    script:
        "scripts/multiple_inputs_one_output.py"

def get_tracking_info_from_subfolderfilename(wildcards):
    subfolder = os.path.split(os.path.split(wildcards.subfolder_filename)[0])[1]
    tracking_info = '3a_tracking_info/' + subfolder + '/output.npy'
    print('this is in the snakemake function 2')
    print(tracking_info)
    return tracking_info


rule test_different_inputs:
    input:
        "1_data/{subfolder_filename}.tif",
        get_tracking_info_from_subfolderfilename,
    output:
        "3b_tracking_images/{subfolder_filename}.tif"
    script:
        "scripts/different_inputs.py"
