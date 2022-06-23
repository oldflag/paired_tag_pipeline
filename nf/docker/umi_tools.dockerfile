# FROM drtools/alpine-conda
FROM continuumio/miniconda3:latest

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="1.0"
LABEL description="This is a custom Docker Image for \
umi_tools, samtools and anndata."

RUN /opt/conda/bin/conda install -c defaults -c conda-forge -c bioconda  --yes --freeze-installed \
    openssl \
    python=3.8 \
    samtools \
    umi_tools \
    anndata \
    awscli \
    && /opt/conda/bin/conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete 
