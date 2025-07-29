# 1
# docker build . -t spectralnanodiamond/btrack --platform=linux/x86_64

# 2
# run docker run -it spectralnanodiamond/btrack
#or
# run docker run -it 124091824019489
# or 
# run docker run -it -p 8888:8888 spectralnanodiamond/btrack

# 3
# docker image push spectralnanodiamond/btrack:latest

# To run jupyter notebook
# run with port option above then
# jupyter notebook --ip 0.0.0.0 --allow-root --no-browser
# then http://127.0.0.1:8888 or http://localhost:8888/

FROM ubuntu:22.04 AS base
# FROM python:3.11-bookworm AS base
WORKDIR /app
ARG MINIFORGE_NAME=Miniforge3
ARG MINIFORGE_VERSION=24.7.1-0
ARG TARGETPLATFORM

ENV CONDA_DIR=/opt/conda
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=${CONDA_DIR}/bin:${PATH}

# 1. Install just enough for conda to work
# 2. Keep $HOME clean (no .wget-hsts file), since HSTS isn't useful in this context
# 3. Install miniforge from GitHub releases
# 4. Apply some cleanup tips from https://jcrist.github.io/conda-docker-tips.html
#    Particularly, we remove pyc and a files. The default install has no js, we can skip that
# 5. Activate base by default when running as any *non-root* user as well
#    Good security practice requires running most workloads as non-root
#    This makes sure any non-root users created also have base activated
#    for their interactive shells.
# 6. Activate base by default when running as root as well
#    The root user is already created, so won't pick up changes to /etc/skel
FROM base AS miniforge-install
RUN apt-get update > /dev/null && \
    apt-get install --no-install-recommends --yes \
        wget bzip2 ca-certificates \
        gcc libsuitesparse-dev libc6-dev liblapack-dev libopenblas-dev liblapack3 \ 
        libopenblas-base libblas-dev libatlas-base-dev pkg-config \ 
        libfontconfig libfontconfig1 qtcreator fontconfig libfontconfig1-dev \ 
        build-essential g++-x86-64-linux-gnu libc6-dev-amd64-cross\
        git \
        tini \
        > /dev/null && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/${MINIFORGE_NAME}-${MINIFORGE_VERSION}-Linux-$(uname -m).sh -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh && \
    conda clean --tarballs --index-cache --packages --yes && \
    find ${CONDA_DIR} -follow -type f -name '*.a' -delete && \
    find ${CONDA_DIR} -follow -type f -name '*.pyc' -delete && \
    conda clean --force-pkgs-dirs --all --yes  && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc
FROM miniforge-install AS checkpoint-cuda


RUN wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb && \
    dpkg -i cuda-keyring_1.1-1_all.deb && \
    apt update

RUN apt-get update && apt-get install -y \
    cuda-cudart-12-3 \
    libcublas-12-3 \
    libcusparse-12-3 \
    libcusolver-12-3 \
    # libcudnn8=8.9.*-1+cuda12.3 \
    && rm -rf /var/lib/apt/lists/*
    # RUN mamba install -y python=3.11 matplotlib jupyterlab natsort jupyter notebook
COPY ./conda_envs_yaml/ /app/conda_envs_yaml/
FROM checkpoint-cuda AS checkpoint-install-base-env
#pyqt napari-arboretum
#pydantic=1.10.14
# COPY /conda_envs_yaml/environment_morphody46_dev.yml environment_morphody46_dev.yml

RUN conda env create -f conda_envs_yaml/environment_morphody46_dev.yml && \
    conda clean --all --yes
RUN conda init bash && . ~/.bashrc && \
    conda activate morphody46
# RUN conda activate morphody46
# echo "conda activate morphodys46" >> ~/.bashrc
# RUN pip install btrack
#[napari]
FROM checkpoint-install-base-env AS checkpoint-cellpose-env


RUN conda env create -y -f conda_envs_yaml/environment_morphody_cellpose2_dev.yml && \
    conda clean --all --yes
RUN conda init bash && . ~/.bashrc

FROM checkpoint-cellpose-env AS checkpoint-btrack-env

RUN conda env create -y -f conda_envs_yaml/environment_morphody_btrack4_dev.yml && \
    conda clean --all --yes
RUN conda init bash && . ~/.bashrc

FROM checkpoint-btrack-env AS final
COPY run_example.smk /app/
CMD [ "/opt/conda/etc/profile.d/conda.sh" ]
