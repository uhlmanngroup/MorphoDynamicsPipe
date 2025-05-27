# MorphoDynamicsPipe
A library that uses snakemake to run a morphodynamics analysis on images of cells or nuclei. 
It includes segmentation via cellpose, StarDist or micro-SAM and tracking via btrack. 
MorpholoDynamics features are calculated with scikit-image, CPDA chord measurement and other custom functions.

## Installation
Installation requires a mamba or conda package manager. Please go here for more details:
https://github.com/conda-forge/miniforge

After installing conda or mamba, clone this repository.
`git clone https://github.com/uhlmanngroup/MorphoDynamicsPipe`

Then navigate to the MorphoDynamicsPipe folder in a terminal or commmand prompt that has conda activated. 

Then run `mamba env create -f conda_envs_yaml/environment_morphody46_dev.yml`
(mamba can also be replaced by conda here)

then

`conda activate morphody46`

Other conda environments (such as for cellpose) will be installed automatically by snakemake in
the default folder for your conda setup. These are currently called
morphody_cellposecpu, morphody_cellpose, morphody_btrackpip, morphody_btrack. 
Not all environments will be installed - only the ones relevant for the hardware you are using 
(e.g. operating system, GPU availability). 

## Data setup
To add data, create a folder called `1_data` inside the main MorphoDynamicsPipe folder. 
Put data in the in the format 
`1_data/folder1/filename1_T=000.tif`, `1_data/folder1/filename1_T=001.tif`, `1_data/folder1/filename1_T=002.tif`, 
`1_data/folder2/filename2_T=000.tif`, `1_data/folder2/filename2_T=001.tif`, `1_data/folder2/filename2_T=002.tif`, 
etc...

where the folder1, folder2, filename1 and filename2 (etc...) can be anything, and filename1 and filename2 (etc...) do not have to be consistent within each subfolder. 
The files must be tiffs with T=[NUMBER].tiff or T=[NUMBER].tif at the end. 
The T= number must run from 0 and increment by one across all the files in each folder. There is no limit on the number of folders or files. 

An example dataset can be found in `MorphoDynamicsPipe/example`.
To run the example data, copy the folder `MorphoDynamicsPipe/example/1_data` to `MorphoDynamicsPipe/1_data`, including the the subfolders and files. 

## Execution
If you are on Windows and do not have a GPU, then open the snakemake file and make a small change. 
Please change `windows_and_cpu_only = False` to `windows_and_cpu_only = True`. 
If you are not on Windows then you can ignore this step. 
If you are on Windows but you have a GPU then you can ignore this step. 

Then run `snakemake -s run_example.smk --cores "all" --sdm conda --keep-going` 
when in the MorphoDynamicsPipe folder.

To change which algorithm is used for segmentation (default: cellpose), open the run_example.smk file, 
comment out all parts of `rule segment_with_cellpose_nucs` and comment in either `rule segment_with_stardist` or `rule segment_with_micro_sam`.

To use nuclei for tracking whole cells, first segment and track nuclei in a separate folder. Then put whole cells in a new MorphoDynamicsPipe folder. 
Comment out `rule segment_with_cellpose_nucs` and comment in `rule segment_with_cellpose_celltracker_with_nucs`. 
Comment out `rule track_with_btrack` and comment in `rule symlink_to_btrack_info` and comment in either the `ln -s` line or the `cp` line (but not both).
Then run as normal. 

If you would like to run the pipeline on another project, simply copy the `.smk` file 
and the `scripts/` folder (including its files) next to another `1_data` 
e.g. `project2/1_data/` , `project2/run_exmaple.smk` and `project2/scripts`. 
Then at the command line, navigate to the project2 folder and run the snakemake command as above.

## Results

Segmentation and tracking results can be viewed in napari (installed with the conda environment) from folder `5_tracking_images_outlines/`.

Morphodynamics results can be found in the folders `4a_instantaneous_cell_morphodynamics/` and `4b_time_averaged_cell_morphodynamics/`.

Example results can be found in the folder `MorphoDynamicsPipe/example/`

macOS seems to produce slightly different segmentation results. 

## Citations
This library relies on: \
Cellpose https://github.com/MouseLand/cellpose \
StarDist https://github.com/stardist/stardist \
Micro-SAM https://github.com/computational-cell-analytics/micro-sam \
btrack https://github.com/quantumjot/btrack \
Scikit-image https://github.com/scikit-image/scikit-image \
CPDA curvature measurement https://ieeexplore.ieee.org/document/4657455 \
among other commonly used packages. 

Example data is taken from Broad Bioimage Benchmark Collection (Broad Institute), 
simulated HL60 cells (from the Cell Tracking Challenge), Accession number BBBC035, Version 1. https://bbbc.broadinstitute.org/BBBC035

This package is part of the PLAST_CELL project. \
More details on this are here: https://plastcell.eu/
