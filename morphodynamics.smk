files = glob_wildcards("1_data/{subfolder_filename}.tif")
both = glob_wildcards("1_data/{subfolder}/{filename}.tif")
both.subfolder
both.filename

rule all:
    input:
        expand("4_regionprops/{subfolder_filename}.csv", subfolder_filename = files.subfolder_filename)

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

rule get_regionprops:
    input:
        "2_segmentation/{subfolder_filename}.tif"
    output:
        "4_regionprops/{subfolder_filename}.csv"
    conda:
        "conda_envs_yaml/environment_StardistSmake_dev.yml"
    shell:
        "python scripts/get_regionprops.py {input} {output}"