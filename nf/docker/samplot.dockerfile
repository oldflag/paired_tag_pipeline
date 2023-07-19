FROM mambaorg/micromamba

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="3.0"
LABEL description="This is a custom Docker Image for samplot"

ADD ./envs/samplot.yaml .
RUN micromamba install -y -n base -f samplot.yaml && \
    micromamba clean --all --yes
RUN micromamba install -c conda-forge procps-ng -n base && \
    micromamba clean --all --yes

# FROM continuumio/miniconda3:latest

# # LABEL 
# LABEL maintainer="hklim@epigenome.us"
# LABEL version="1.0"
# LABEL description="This is custom Docker Image for \
# samplot"

# RUN /opt/conda/bin/conda install -c defaults -c conda-forge -c bioconda  --yes --freeze-installed \
#     samplot \
#     samtools \
#     bedtools \
#     tabix \
#     awscli \
#     && /opt/conda/bin/conda clean -afy \
#     && find /opt/conda/ -follow -type f -name '*.a' -delete \
#     && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
#     && find /opt/conda/ -follow -type f -name '*.js.map' -delete 
