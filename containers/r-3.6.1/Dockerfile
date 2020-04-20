FROM nfcore/base:1.9
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing statistics requirements for the nfcore/samba pipeline"

COPY containers/r-3.6.1/environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-r-3.6.1/bin:$PATH
