FROM mambaorg/micromamba

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="3.0"
LABEL description="This is a custom Docker Image for samtools, pysam, etc."

ADD ./envs/bamtofrag.yaml .
RUN micromamba install -y -n base -f bamtofrag.yaml && \
    micromamba clean --all --yes
