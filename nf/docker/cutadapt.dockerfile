# FROM drtools/alpine-conda
FROM continuumio/miniconda3:latest

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="1.0"
LABEL description="This is a custom Docker Image for \
cutadapt."

RUN /opt/conda/bin/conda install -c defaults -c conda-forge -c bioconda --yes --freeze-installed \
    cutadapt \
    nomkl \
    awscli \
    && /opt/conda/bin/conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete 
