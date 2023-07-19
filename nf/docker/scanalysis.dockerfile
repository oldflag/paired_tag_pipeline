FROM mambaorg/micromamba

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="3.0"
LABEL description="This is a custom Docker Image for scanpy, seabone, anndata."

ADD ./envs/scanalysis.yaml .
RUN micromamba install -y -n base -f scanalysis.yaml && \
    micromamba clean --all --yes
RUN micromamba install -c conda-forge procps-ng -n base && \
    micromamba clean --all --yes
# FROM drtools/alpine-conda
# FROM continuumio/miniconda3:latest

# # LABEL 
# LABEL maintainer="hklim@epigenome.us"
# LABEL version="1.0"
# LABEL description="This is a custom Docker Image for \
# scanpy, seabone, anndata."

# RUN /opt/conda/bin/conda install -c defaults -c conda-forge -c bioconda --yes --freeze-installed \
#     seaborn \
#     anndata==0.8.0 \
#     scanpy \
#     pandas \
#     numpy==1.21 \
#     leidenalg \
#     harmonypy \
#     nomkl \
#     awscli \
#     && /opt/conda/bin/conda clean -afy \
#     && find /opt/conda/ -follow -type f -name '*.a' -delete \
#     && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
#     && find /opt/conda/ -follow -type f -name '*.js.map' -delete 
