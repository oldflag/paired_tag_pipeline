# FROM drtools/alpine-conda
FROM continuumio/miniconda3:latest

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="1.5"
LABEL description="scikit-bio, seaborn, samtools,etc for paired-tag process"

RUN /opt/conda/bin/conda install -c defaults -c conda-forge -c bioconda --yes --freeze-installed \
    scikit-bio \
    seaborn \
    samtools \
    pysam \
    distance \
    logomaker \
    awscli \
    && /opt/conda/bin/conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete 