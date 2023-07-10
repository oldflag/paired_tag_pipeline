FROM mambaorg/micromamba

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="3.0"
LABEL description="This is a custom Docker Image for \
cutadapt"

ADD ./envs/cutadapt.yaml .
RUN micromamba install -y -n base -f cutadapt.yaml && \
    micromamba clean --all --yes
