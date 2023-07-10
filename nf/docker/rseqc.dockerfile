FROM mambaorg/micromamba

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="3.0"
LABEL description="This is a custom Docker Image for rna-seqc"

ADD ./envs/rseqc.yaml .
RUN micromamba install -y -n base -f rseqc.yaml && \
    micromamba clean --all --yes

# FROM drtools/alpine-conda
# FROM continuumio/miniconda3:latest

# # LABEL 
# LABEL maintainer="hklim@epigenome.us"
# LABEL version="1.0"
# LABEL description="This is a custom Docker Image for \
# rna-seqc."

# RUN /opt/conda/bin/conda install -c conda-forge -c bioconda --yes --freeze-installed \
#     rna-seqc \
#     nomkl \
#     awscli \
#     && /opt/conda/bin/conda clean -afy \
#     && find /opt/conda/ -follow -type f -name '*.a' -delete \
#     && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
#     && find /opt/conda/ -follow -type f -name '*.js.map' -delete 
