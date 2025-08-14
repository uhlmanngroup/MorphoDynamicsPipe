#!/bin/bash
set -e

echo "[INFO] Checking for NVIDIA GPU availability..."

if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
  echo "[INFO] GPU found. Launching 'app-gpu' service."
  docker compose up morphodys46 "$@"
else
  echo "[WARNING] No GPU detected. Launching 'app-cpu' service."
  docker compose up morphodys46-cpu "$@"
fi