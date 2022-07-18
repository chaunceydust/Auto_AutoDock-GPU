FROM nvidia/cuda:11.3.1-devel-ubuntu20.04
MAINTAINER Jongseo_Park jongseopark@gm.gist.ac.kr

# Initial setup
RUN rm /etc/apt/sources.list.d/cuda.list \
    && apt-key del 7fa2af80

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Seoul
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone \
    && apt-get update \
    && apt-get install -y \
    wget \
    git \
    python3 \
    python3-pip \
    unzip \
    cmake \
    tzdata \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Conda setup
RUN wget -q -P /tmp https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniconda3-latest-Linux-x86_64.sh \
    && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
    && echo 'export PATH="/opt/conda/bin/":$PATH' >> ~/.bashrc \
    && /bin/bash -c "source ~/.bashrc"


# Git clone (auto_autodock-gpu / openbabel / autodock_gpu)
WORKDIR /opt/
RUN git clone https://github.com/jongseo-park/Auto_AutoDock-GPU \
    && git clone https://github.com/ccsb-scripps/AutoDock-GPU \
    && wget https://github.com/openbabel/openbabel/archive/refs/tags/openbabel-2-4-1.tar.gz \
    && tar -zxvf openbabel-2-4-1.tar.gz \
    && rm -rf openbabel-2-4-1.tar.gz

# autodock-gpu setup
WORKDIR /opt/AutoDock-GPU
RUN export GPU_INCLUDE_PATH="/usr/local/cuda/include/" \
    && export GPU_LIBRARY_PATH="/usr/local/cuda/lib64" \ 
    && make DEVICE=CUDA NUMWI=128

# openbabel setup
WORKDIR /opt/openbabel-openbabel-2-4-1
RUN mkdir build

WORKDIR /opt/openbabel-openbabel-2-4-1/build 
RUN cmake .. -DCMAKE_INSTALL_PREFIX=/opt/openbabel/ \
    && make -j 4 \
    && make install


# Generate conda envs
WORKDIR /opt/Auto_AutoDock-GPU
ENV PATH="/opt/conda/bin:$PATH"
RUN conda env create -f requirements.yml \ 
    && conda clean --all


# Finalize
WORKDIR /home/run/
RUN echo 'export PATH="/opt/openbabel/bin/":$PATH' >> ~/.bashrc \
    && echo "conda activate autodock_gpu" >> ~/.bashrc \