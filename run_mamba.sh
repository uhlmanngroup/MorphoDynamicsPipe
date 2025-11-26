#!/bin/bash
#chmod +x run_mamba.sh
#source run_mamba.sh
set -e


if command -v mamba >/dev/null 2>&1; then
  echo "[INFO] 'mamba' found."
  PKG_MGR="mamba"
  # Initialize shell for mamba via hook only (no fallback to conda base)
  hook_cmd="$(mamba shell hook --shell bash 2>/dev/null || true)"
  SHELL_INITIALIZED=false
  if [ -n "$hook_cmd" ]; then
    eval "$hook_cmd"
    SHELL_INITIALIZED=true
  else
    echo "[WARNING] mamba shell hook not available in this shell. Will avoid 'mamba activate' and use 'mamba run' for in-env commands."
  fi
  # check if env exists before creating (use mamba env list)
  if mamba env list 2>/dev/null | awk '{print $1}' | grep -qx 'morphody50'; then
    echo "[INFO] Environment 'morphody50' already exists."
    if [ "$SHELL_INITIALIZED" = true ]; then
      echo "[INFO] Activating 'morphody50'..."
      mamba activate morphody50
      snakemake -s run_example.smk --cores 4 --keep-going
    else
      echo "[INFO] Shell not initialized; skipping activation. Use 'mamba run -n morphody50 <cmd>' to run commands inside the env."
    fi
  else
    echo "[INFO] Creating environment 'morphody50'..."
    mamba create -c conda-forge -c nodefaults -n morphody50 python==3.11.11 pip==25.3 -y
    if [ "$SHELL_INITIALIZED" = true ]; then
      mamba activate morphody50
      pip install -r requirements-base.txt
      snakemake -s run_example.smk --cores 4 --keep-going
    else
      echo "[INFO] Installing requirements into 'morphody50' using 'mamba run'..."
#        mamba run -n morphody50 pip install -r requirements-base.txt
    fi

  fi
else
  echo "[ERROR] 'mamba' not found. Please install Mamba to proceed."
  exit 1
fi