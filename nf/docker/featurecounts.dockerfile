# FROM drtools/alpine-conda
FROM continuumio/miniconda3:latest

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="1.0"
LABEL description="This is a custom Docker Image for \
featurecount, samtools and pysam."

RUN /opt/conda/bin/conda install  -c conda-forge -c bioconda -c defaults --yes --freeze-installed \
    subread \
    openssl \
    samtools=1.12 \
    pysam \
    nomkl \
    awscli \
    && /opt/conda/bin/conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete 
