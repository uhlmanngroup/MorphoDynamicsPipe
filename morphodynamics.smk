rule segment_with_stardist:
    input:
        "../data/301432316.tif"
    output:
        "results/1.npy"
    conda:
        "conda_envs_yaml/environment_StardistSmake_dev.yml"
    script:
        "scripts/change_format.py"
 #   shell:
 #       "python scripts/change_format.py {input} {output}"



#    shell: "cp {input} {output}"