FROM mambaorg/micromamba

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="3.0"
LABEL description="This is a custom Docker Image for scanpy, seabone, anndata."

ADD ./envs/skbio.yaml .
RUN micromamba install -y -n base -f skbio.yaml && \
    micromamba clean --all --yes
RUN micromamba install -c conda-forge procps-ng -n base && \
    micromamba clean --all --yes
RUN micromamba install -c conda-forge ncurses -n base && \
    micromamba clean --all --yes

    
# && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# # FROM drtools/alpine-conda
# FROM continuumio/miniconda3:latest

# # LABEL 
# LABEL maintainer="hklim@epigenome.us"
# LABEL version="1.5"
# LABEL description="scikit-bio, seaborn, samtools,etc for paired-tag process"

# RUN /opt/conda/bin/conda install -c defaults -c conda-forge -c bioconda --yes --freeze-installed \
#     scikit-bio \
#     seaborn \
#     samtools \
#     pysam \
#     distance \
#     logomaker \
#     awscli \
#     && /opt/conda/bin/conda clean -afy \
#     && find /opt/conda/ -follow -type f -name '*.a' -delete \
#     && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
#     && find /opt/conda/ -follow -type f -name '*.js.map' -delete 