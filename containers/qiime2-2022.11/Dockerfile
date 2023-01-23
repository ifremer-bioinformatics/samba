FROM --platform=linux/x86_64 mambaorg/micromamba
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing QIIME2 v2022.11 plugin for QIIME2 requirements for the samba pipeline"

COPY qiime2-2022.11.yml /
RUN micromamba create -f /qiime2-2022.11.yml

SHELL ["micromamba", "run", "-n", "qiime2-2022.11", "/bin/bash", "-c"]
RUN mkdir -p /opt
WORKDIR /opt
RUN git clone https://github.com/cduvallet/q2-dbotu
WORKDIR /opt/q2-dbotu
USER root
RUN micromamba install -c cduvallet -c conda-forge procps-ng biopython dbotu && python setup.py install
RUN R --vanilla -e 'install.packages("devtools", repos="https://mirror.ibcp.fr/pub/CRAN/")'
RUN R --vanilla -e 'devtools::install_github("donaldtmcknight/microDecon")'
RUN micromamba install -c bioconda seqtk
RUN micromamba clean -a
ENV PATH /opt/conda/envs/qiime2-2022.11/bin:$PATH
