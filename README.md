# MorphoDynamicsPipe
A library to use snakemake to run a morphodynamics analysis on cells

run `mamba env create -f conda_envs_yaml/environment_morphody8_dev.yml`

`conda activate morphody8`
Then run `snakemake -s morphodynamics.smk --cores='all' --sdm conda --conda-frontend mamba --keep-going` 
when in the folder


To add data, create a folder called `1_data` inside the main MorphoDynamicsPipe folder. 
Put data in the in the format `1_data/folder1/file1.tif`, `1_data/folder1/file2.tif`