files = glob_wildcards("../data/{code}.tif")

rule all:
    input:
        expand("results/{code}.npy", code = files.code)


rule segment_with_stardist:
    input:
        "../data/{code}.tif"
    output:
        "results/{code}.npy"
    conda:
        "conda_envs_yaml/environment_StardistSmake_dev.yml"
    shell:
        "python scripts/run_stardist.py {input} {output}"