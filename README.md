# MorphoDynamicsPipe
A library to use snakemake to run a morphodynamics analysis on cells

run mamba env create -f environment_Smake_dev.yml
conda activate Smake0
Then run snakemake -s morphodynamics.smk --cores='all' when in the folder
