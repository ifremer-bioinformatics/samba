FROM --platform=linux/x86_64 mambaorg/micromamba
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing PICRUSt2 v2.5.1 for the samba pipeline"

COPY picrust2-2.5.1.yml /
RUN micromamba create -f /picrust2-2.5.1.yml

USER root
SHELL ["micromamba", "run", "-n", "picrust2-2.5.1", "/bin/bash", "-c"]
RUN mkdir -p /opt
WORKDIR /opt
RUN wget https://github.com/picrust/picrust2/archive/v2.5.1.tar.gz
RUN tar -xvzf  v2.5.1.tar.gz
WORKDIR /opt/picrust2-2.5.1/
RUN pip install --editable .

RUN micromamba clean -a
ENV PATH /opt/conda/envs/picrust2-2.5.1/bin:$PATH
