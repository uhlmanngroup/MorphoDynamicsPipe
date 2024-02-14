rule segment_with_stardist:
    input:
        "../data/301432316.tif"
    output:
        "results/1.npy"
    conda:
        "conda_envs_yaml/environment_StardistSmake_dev.yml"
    shell:
        "python scripts/change_format.py {input} {output}"