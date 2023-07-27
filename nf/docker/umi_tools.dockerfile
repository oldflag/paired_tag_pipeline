FROM mambaorg/micromamba

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="3.0"
LABEL description="This is a custom Docker Image for umi_tools, samtools and anndata."

ADD ./envs/umi_tools.yaml .
RUN micromamba install -y -n base -f umi_tools.yaml && \
    micromamba clean --all --yes
RUN micromamba install -c conda-forge procps-ng -n base && \
    micromamba clean --all --yes
RUN micromamba install -c conda-forge ncurses -n base && \
    micromamba clean --all --yes


# FROM drtools/alpine-conda
# FROM continuumio/miniconda3:latest

# # LABEL 
# LABEL maintainer="hklim@epigenome.us"
# LABEL version="1.0"
# LABEL description="This is a custom Docker Image for \
# umi_tools, samtools and anndata."

# RUN /opt/conda/bin/conda install -c defaults -c conda-forge -c bioconda  --yes --freeze-installed \
#     openssl \
#     python=3.8 \
#     samtools \
#     umi_tools \
#     anndata \
#     awscli \
#     && /opt/conda/bin/conda clean -afy \
#     && find /opt/conda/ -follow -type f -name '*.a' -delete \
#     && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
#     && find /opt/conda/ -follow -type f -name '*.js.map' -delete 
