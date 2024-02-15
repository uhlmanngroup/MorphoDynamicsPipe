files = glob_wildcards("../data/{subfolder_filename}.tif")

rule all:
    input:
        expand("results/{subfolder_filename}.tif", subfolder_filename = files.subfolder_filename)


rule segment_with_stardist:
    input:
        "../data/{subfolder_filename}.tif"
    output:
        "results/{subfolder_filename}.tif"
    conda:
        "conda_envs_yaml/environment_StardistSmake_dev.yml"
    shell:
        "python scripts/run_stardist.py {input} {output}"