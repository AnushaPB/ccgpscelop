# Use rocker/geospatial as the base image
FROM rocker/geospatial:latest

# Install system dependencies for Python and FEEMS
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-dev \
    libgeos-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda to manage Python packages
ENV MINICONDA_VERSION 4.7.12.1
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_${MINICONDA_VERSION}-Linux-x86_64.sh -O ~/miniconda.sh \
    && /bin/bash ~/miniconda.sh -b -p /opt/conda \
    && rm ~/miniconda.sh \
    && /opt/conda/bin/conda clean -tipsy

# Add Conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Create a Conda environment and install FEEMS
RUN conda create -n feems_e python=3.8.3 -y
SHELL ["conda", "run", "-n", "feems_e", "/bin/bash", "-c"]

# Activating the environment and installing FEEMS and its dependencies
RUN conda install -c bioconda feems -c conda-forge -y