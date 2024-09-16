# run docker build . -t spectralnanodiamond/btrack --platform=linux/x86_64
# run docker run -it 124091824019489
# docker image push spectralnanodiamond/btrack:latest

#FROM ubuntu:focal-20240530
FROM ubuntu:22.04

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
RUN mamba install -y python=3.11.9 pydantic=1.10.14 btrack

ENTRYPOINT ["tini", "--"]
CMD [ "/bin/bash" ]
