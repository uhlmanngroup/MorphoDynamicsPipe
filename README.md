# MorphoDynamicsPipe
A library to use snakemake to run a morphodynamics analysis on cells.

## Installation
run `mamba env create -f conda_envs_yaml/environment_morphody35_dev.yml`
`conda activate morphody35`

## Data setup
To add data, create a folder called `1_data` inside the main MorphoDynamicsPipe folder. 
Put data in the in the format 
`1_data/folder1/file1_T=0.tif`, `1_data/folder1/file2_T=1.tif`, `1_data/folder1/file3_T=2.tif`, 
`1_data/folder2/file1_T=0.tif`, `1_data/folder2/file2_T=1.tif`, `1_data/folder2/file3_T=2.tif`, 

where the folder and file names can be anything, but the files must be tiffs with T=[NUMBER].tiff or T=[NUMBER].tif at the end. 
The T= number must run from 0 and increment by one across all the files in each folder. 

## Execution
Then run `snakemake -s run_example.smk --cores='all' --sdm conda --conda-frontend mamba --keep-going` 
when in the MorphoDynamicsPipe folder.
