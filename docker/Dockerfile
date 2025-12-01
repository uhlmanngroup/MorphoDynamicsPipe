#FROM nvidia/cuda:12.1.1-runtime-ubuntu22.04 AS base
#FROM pytorch/pytorch:2.4.1-cuda12.1-cudnn9-devel AS base
FROM pytorch/pytorch:2.8.0-cuda12.6-cudnn9-devel AS base
WORKDIR /app
COPY ./requirements-base.txt /app/requirements.txt
RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends \
    git curl ca-certificates build-essential graphviz \
    libgl1-mesa-dri libegl1-mesa libgles2-mesa xvfb \
    libsm6 libxext6 libxrender1 ffmpeg \
    && rm -rf /var/lib/apt/lists/*
    #python3 python3-pip python3-venv python3-dev
#RUN wget https://www.python.org/ftp/python/3.11.11/Python-3.11.11.tgz \
#    && tar xvf Python-3.11.11.tgz \
#    && cd Python-3.11.11 \
#    && ./configure --enable-optimizations \
#    && make -j$(nproc) \
#    && make install \
#    && cd .. \
#    && rm -rf Python-3.11.11 Python-3.11.11.tgz
FROM base AS venv_install
#RUN apt-get update && apt-get upgrade -y && apt-get install -y python3-venv && rm -rf /var/lib/apt/lists/*
#RUN pip install virtualenv
#RUN python -m ensurepip --upgrade
#RUN python -m venv venv && . ./venv/bin/activate && pip install --no-cache-dir -r ./requirements.txt
RUN conda init
RUN python --version
RUN pip install --no-cache-dir -r ./requirements.txt
FROM venv_install AS final
RUN apt clean && rm -rf /var/lib/apt/lists/*

