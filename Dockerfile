FROM python:3.11-slim AS base

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# 1. Installer les dépendances système de base
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    gnupg \
    ca-certificates \
    libgl1 \
    libglib2.0-0 \
    git \
    pkg-config \
    libx11-dev \
    libxi6 \
    libxrender1 \
    libxtst6 \
    libnss3 \
    libatk1.0-0 \
    libatk-bridge2.0-0 \
    libgtk-3-0 \
    && rm -rf /var/lib/apt/lists/*

FROM base AS builder
# 2. Ajouter le dépôt CUDA officiel de NVIDIA (runtime only, pas toolkit complet)
RUN mkdir -p /etc/apt/keyrings && \
    curl -fsSL https://developer.download.nvidia.com/compute/cuda/repos/debian12/x86_64/3bf863cc.pub \
    | gpg --dearmor -o /etc/apt/keyrings/nvidia.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/nvidia.gpg] https://developer.download.nvidia.com/compute/cuda/repos/debian12/x86_64 /" \
    > /etc/apt/sources.list.d/cuda.list

FROM builder AS checkpoint-cuda
RUN apt-get update && apt-get install -y \
    cuda-cudart-12-3 \
    libcublas-12-3 \
    libcusparse-12-3 \
    libcusolver-12-3 \
    # libcudnn8=8.9.*-1+cuda12.3 \
    && rm -rf /var/lib/apt/lists/*

# 3. Définir les variables d'environnement CUDA
FROM checkpoint-cuda AS checkpoint-install-base
ENV PATH=/usr/local/cuda/bin:${PATH}
# ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64:/usr/local/cuda/lib:${LD_LIBRARY_PATH}

# 4. Installer pip + requirements
WORKDIR /app
COPY requirements-base.txt .
RUN pip install --upgrade pip && pip install -r requirements-base.txt

FROM checkpoint-install-base AS checkpoint-install-dev
COPY requirements-dev.txt .
RUN pip install --upgrade pip && pip install -r requirements-dev.txt

FROM checkpoint-install-dev AS checkpoint-install-extended
COPY requirements-extended.txt .
RUN pip install --upgrade pip && pip install -r requirements-extended.txt

# FROM passing AS notpassing
# COPY requirements_not_passing.txt .
# RUN pip install --upgrade pip && pip install -r requirements_not_passing.txt

FROM checkpoint-install-extended AS verification
RUN pip check

FROM verification AS final
# 5. Définir la commande par défaut
CMD ["python3"]
