# FROM drtools/alpine-conda
FROM continuumio/miniconda3:latest

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="1.0"
LABEL description="This is a custom Docker Image for \
scanpy, seabone, anndata."

RUN /opt/conda/bin/conda install -c defaults -c conda-forge -c bioconda --yes --freeze-installed \
    seaborn \
    anndata==0.8.0 \
    scanpy \
    pandas \
    numpy \
    nomkl \
    awscli \
    && /opt/conda/bin/conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete 
