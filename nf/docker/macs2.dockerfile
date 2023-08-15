FROM mambaorg/micromamba

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="3.0"
LABEL description="This is a custom Docker Image for \
macs2, samtools and bedtools."
USER root
RUN apt-get update && apt-get install sudo
RUN apt-get install -y git
RUN apt-get install -y python3-pip
USER mambauser
ADD ./envs/macs2.yaml .
RUN micromamba install -y -n base -f macs2.yaml && \
    micromamba clean --all --yes
RUN micromamba install -c conda-forge procps-ng -n base && \
    micromamba clean --all --yes
RUN micromamba install -c conda-forge ncurses -n base && \
    micromamba clean --all --yes
RUN git clone --recurse-submodules https://gitlab.com/epigenomeus_public/MACS.git ~/MACS
# RUN pip install --upgrade pip
# RUN pip install --upgrade setuptools
RUN DEB_PYTHON_INSTALL_LAYOUT=deb_system pip install ~/MACS
RUN alias macs2=macs3

