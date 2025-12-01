#!/bin/bash
#docker login
#docker buildx create --use
#docker buildx build --platform linux/amd64 -t spectralnanodiamond/morphodynamicspipe:latest --push .
#docker build --platform linux/amd64 -t spectralnanodiamond/morphodynamicspipe -f docker/Dockerfile .
#chmod +x run_docker.sh
#source run_docker.sh
#python -c "import torch; print(torch.__version__); print('CUDA:', torch.version.cuda); print('cuda available:', torch.cuda.is_available())"

set -e

# Check that docker CLI is available
if ! command -v docker >/dev/null 2>&1; then
  echo "[ERROR] Docker CLI not found. Please install Docker and ensure 'docker' is installed and on your PATH."
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
