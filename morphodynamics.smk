rule segment_with_stardist:
    input:
        "../data/301432316.tif"
    output:
        "results/301432316.tif"
    shell: "cp {input} {output}"