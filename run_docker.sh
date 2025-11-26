#!/bin/bash
#chmod +x run_docker.sh
#source run_docker.sh
set -e

# Check that docker CLI is available
if ! command -v docker >/dev/null 2>&1; then
  echo "[ERROR] Docker CLI not found. Please install Docker and ensure 'docker' is on your PATH."
  echo "[HINT] Installation: https://docs.docker.com/get-docker/"
  exit 1
fi

if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
  echo "[INFO] GPU found. Launching 'app-gpu' service."
  docker compose up morphody50 "$@"
else
  echo "[WARNING] No GPU detected. Launching 'app-cpu' service."
  docker compose up morphody50-cpu "$@"
fi
