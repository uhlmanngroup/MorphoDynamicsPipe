#I want to create a shell script that does the same as the docker compose up commands and the dockerfile, but is a singularity command
#!/bin/bash
#chmod +x run_singularity.sh
#source run_singularity.sh
set -e

singularity run -i docker://spectralnanodiamond/morphodynamicspipe:latest sh
#singularity pull docker://spectralnanodiamond/morphodynamicspipe:latest
