# MorphoDynamicsPipe
A library to use snakemake to run a morphodynamics analysis on cells

run `mamba env create -f conda_envs_yaml/environment_morphody33_dev.yml`

`conda activate morphody33`
Then run `snakemake -s morphodynamics.smk --cores='all' --sdm conda --conda-frontend mamba --keep-going` 
when in the folder


To add data, create a folder called `1_data` inside the main MorphoDynamicsPipe folder. 
Put data in the in the format 
`1_data/folder1/file1_T=0.tif`, `1_data/folder1/file2_T=1.tif`, `1_data/folder1/file3_T=2.tif`, 
`1_data/folder2/file1_T=0.tif`, `1_data/folder2/file2_T=1.tif`, `1_data/folder2/file3_T=2.tif`, 

where the folder and file names can be anything, but the files must be tiffs with T=[NUMBER].tiff at the end. 
The T= number must run from 0 and increment by one across all the files in each folder. 
