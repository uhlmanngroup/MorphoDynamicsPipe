files = glob_wildcards("1_data/{subfolder_filename}.tif")

rule all:
    input:
        expand("2_segmentation/{subfolder_filename}.tif", subfolder_filename = files.subfolder_filename)

rule segment_with_stardist:
    input:
        "1_data/{subfolder_filename}.tif"
    output:
        "2_segmentation/{subfolder_filename}.tif"
    retries: 10
    conda:
        "conda_envs_yaml/environment_StardistSmake_dev.yml"
    shell:
        "python scripts/run_stardist.py {input} {output}"

#rule track_with_btrack:
#    input:
#        "../data/{subfolder_filename}.tif"
#    output:
#        "results/{subfolder_filename}.tif"
#    conda:
#        "conda_envs_yaml/environment_BtrackSmake_dev.yml"
#    shell:
#        "python scripts/run_btrack.py {input} {output}"
#```

rule get_regionprops:
    input:
        "2_segmentation/{subfolder_filename}.tif"
    output:
        "4_regionprops/{subfolder_filename}.csv"
    conda:
        "conda_envs_yaml/environment_StardistSmake_dev.yml"
    shell:
        "python scripts/get_regionprops.py {input} {output}"