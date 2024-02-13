rule segment_with_stardist:
    input:
        "../data/301432316.tif"
    output:
        "results/1.npy"
    shell:
        "python scripts/change_format.py {input} {output}"



#    shell: "cp {input} {output}"