FROM python:3.11-slim

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

# 2. Ajouter le dépôt CUDA officiel de NVIDIA (runtime only, pas toolkit complet)
RUN mkdir -p /etc/apt/keyrings && \
    curl -fsSL https://developer.download.nvidia.com/compute/cuda/repos/debian12/x86_64/3bf863cc.pub \
    | gpg --dearmor -o /etc/apt/keyrings/nvidia.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/nvidia.gpg] https://developer.download.nvidia.com/compute/cuda/repos/debian12/x86_64 /" \
    > /etc/apt/sources.list.d/cuda.list && \
    apt-get update && apt-get install -y \
    cuda-cudart-12-3 \
    libcublas-12-3 \
    libcusparse-12-3 \
    libcusolver-12-3 \
    libcudnn8=8.9.*-1+cuda12.3 \
    && rm -rf /var/lib/apt/lists/*

# 3. Définir les variables d'environnement CUDA
ENV PATH=/usr/local/cuda/bin:${PATH}
ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64:/usr/local/cuda/lib:${LD_LIBRARY_PATH}

# 4. Installer pip + requirements
WORKDIR /app
COPY requirements.txt .
RUN pip install --upgrade pip && pip install -r requirements.txt

# 5. Définir la commande par défaut
CMD ["python3"]
