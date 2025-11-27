#!/bin/bash
#chmod +x run_conda.sh
#source run_conda.sh
set -e

#echo "[INFO] Detecting OS and architecture..."
#echo "[INFO] Checking for NVIDIA GPU availability..."
#os="$(uname -s)"
#arch="$(uname -m)"
   
if command -v conda >/dev/null 2>&1; then
  echo "[INFO] 'conda' found."
  PKG_MGR="conda"
  # Initialize shell for conda if possible
  hook_cmd="$(conda shell.bash hook 2>/dev/null || true)"
  if [ -n "$hook_cmd" ]; then
    eval "$hook_cmd"
  else
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    if [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
      # shellcheck source=/dev/null
      source "$CONDA_BASE/etc/profile.d/conda.sh"
    fi
  fi
  if conda info --envs 2>/dev/null | awk '{print $1}' | grep -qx 'morphody50'; then
    echo "[INFO] Environment 'morphody50' already exists. Activating..."
    conda activate morphody50
    snakemake -s run_example.smk --cores 4 --keep-going
  else
    echo "[INFO] Creating environment 'morphody50'..."
    conda create -c conda-forge -c nodefaults -n morphody50 python==3.11.13 pip==25.1.1 -y
    conda activate morphody50
    pip install -r requirements-base.txt
    snakemake -s run_example.smk --cores 4 --keep-going
  fi
else
  echo "[WARNING] 'conda' not found."
  PKG_MGR=""
fi