FROM mambaorg/micromamba

# LABEL 
LABEL maintainer="hklim@epigenome.us"
LABEL version="3.0"
LABEL description="This is a custom Docker Image for \
macs2, samtools and bedtools."

USER root
RUN apt-get update && apt-get install sudo
RUN apt-get install -y git
RUN apt-get install -y python3
RUN apt-get install -y python3-pip
RUN apt-get install python3-setuptools
USER $MAMBA_USER
ADD ./envs/macs2.yaml .
RUN micromamba install -y -n base -f macs2.yaml && \
    micromamba clean --all --yes
RUN micromamba install -c conda-forge procps-ng -n base && \
    micromamba clean --all --yes
RUN micromamba install -c conda-forge ncurses -n base && \
    micromamba clean --all --yes
RUN python3
RUN pip install cython
RUN pip install numpy==1.22
RUN git clone --recurse-submodules https://gitlab.com/epigenomeus_public/MACS.git ~/MACS
# RUN pip install --upgrade pip
# RUN pip install --upgrade --use-pep517 setuptools
RUN pip install git+https://github.com/huggingface/transformers
RUN pip install cykhash>=2.0
RUN pip install hmmlearn>=0.3
RUN pip install ~/MACS
#RUN cd ~/MACS && python3 setup.py install
#RUN echo "alias macs2='macs3'" >>~/.bashrc
RUN ln -s ~/.local/bin/macs3 ~/.local/bin/macs2
ENV PATH="$PATH:~/.local/bin"
#CMD /bin/bash -c "source ~/.bashrc"
#RUN alias macs2='macs3'

