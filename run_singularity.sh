#I want to create a shell script that does the same as the docker compose up commands and the dockerfile, but is a singularity command
#singularity pull morphodynamicspipe.sif docker://spectralnanodiamond/morphodynamicspipe:latest
#salloc -t 8:00:00 --mem=128G -c 16 --gres=gpu:1
#sbatch -t 24:00:00 --mem=64G -c 16 --gres=gpu:1 --wrap="bash run_singularity.sh"
#add run.smk to bind mounts and change snakemake command to use run.smk if necessary
#singularity run -B ./1_data:/app/1_data   -B ./2_segmentation:/app/2_segmentation   -B ./3a_tracking_info:/app/3a_tracking_info   -B ./3b_tracking_images:/app/3b_tracking_images   -B ./3c_tracking_images_filtered:/app/3c_tracking_images_filtered   -B ./4a_instantaneous_cell_morphodynamics:/app/4a_instantaneous_cell_morphodynamics   -B ./4b_time_averaged_cell_morphodynamics:/app/4b_time_averaged_cell_morphodynamics   -B ./5_tracking_images_outlines:/app/5_tracking_images_outlines   -B ./scripts:/app/scripts   -B ./run_example.smk:/app/run_example.smk     -B ./run.smk:/app/run.smk   -B ./btrack_cell_config.json:/app/btrack_cell_config.json   -B ./.snakemake:/app/.snakemake   -B ./maximum_common_time.txt:/app/maximum_common_time.txt /nfs/research/uhlmann/bwoodhams/git/singularity/morphodynamicspipe.sif
#--rerun-incomplete can solve some issues on the singularity runs where snakemake doesn't complete all steps in one go

#!/bin/bash
# chmod +x run_singularity.sh
# ./run_singularity.sh
set -e

CACHE_DIR="${SINGULARITY_CACHEDIR:-$HOME/.singularity/cache}"
SIF_FILE="$CACHE_DIR/morphodynamicspipe.sif"
IMAGE_URI="docker://spectralnanodiamond/morphodynamicspipe:latest"

# Check singularity is available
if ! command -v singularity >/dev/null 2>&1; then
  echo "[ERROR] Singularity not found. Please install Singularity."
  exit 1
fi

# Create cache directory if it doesn't exist
mkdir -p "$CACHE_DIR"

# Pull image if .sif doesn't exist or if --update flag is provided
if [ ! -f "$SIF_FILE" ] || [[ "$*" == *"--update"* ]]; then
  echo "[INFO] Pulling $IMAGE_URI to $SIF_FILE..."
  singularity pull --arch amd64 "$SIF_FILE" "$IMAGE_URI"
else
  echo "[INFO] Using cached image: $SIF_FILE"
fi

# Detect GPU and set flags
GPU_FLAG=""
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
  echo "[INFO] GPU found. Enabling --nv flag."
  GPU_FLAG="--nv"
else
  echo "[WARNING] No GPU detected. Running CPU mode."
fi

# Create required directories if they don't exist
DIRS=("2_segmentation" "3a_tracking_info" "3b_tracking_images" "3c_tracking_images_filtered" "4a_instantaneous_cell_morphodynamics" "4b_time_averaged_cell_morphodynamics" "5_tracking_images_outlines" ".snakemake")
for dir in "${DIRS[@]}"; do
  if [ ! -d "$dir" ]; then
    echo "[INFO] Creating directory: $dir"
    mkdir -p "$dir"
  fi
done

# Create maximum_common_time.txt if it doesn't exist
if [ ! -f "maximum_common_time.txt" ]; then
  echo "[INFO] Creating file: maximum_common_time.txt"
  echo "1000000000000000000000000000000" > maximum_common_time.txt
fi

# Run singularity with bind mounts (equivalent to docker volumes)
echo "[INFO] Starting singularity container..."
singularity run $GPU_FLAG -B ./1_data:/app/1_data \
  -B ./2_segmentation:/app/2_segmentation \
  -B ./3a_tracking_info:/app/3a_tracking_info \
  -B ./3b_tracking_images:/app/3b_tracking_images \
  -B ./3c_tracking_images_filtered:/app/3c_tracking_images_filtered \
  -B ./4a_instantaneous_cell_morphodynamics:/app/4a_instantaneous_cell_morphodynamics \
  -B ./4b_time_averaged_cell_morphodynamics:/app/4b_time_averaged_cell_morphodynamics \
  -B ./5_tracking_images_outlines:/app/5_tracking_images_outlines \
  -B ./scripts:/app/scripts \
  -B ./run_example.smk:/app/run_example.smk \
  -B ./btrack_cell_config.json:/app/btrack_cell_config.json \
  -B ./.snakemake:/app/.snakemake \
  -B ./maximum_common_time.txt:/app/maximum_common_time.txt \
  "$SIF_FILE" bash -c "cd /app && pwd && ls -l && snakemake -s run_example.smk --cores 4 --keep-going"

