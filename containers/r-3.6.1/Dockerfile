FROM nfcore/base:1.9
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing statistics requirements for the samba pipeline"

RUN apt-get update
RUN apt-get install -y librsvg2-dev
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/r-3.6.1/bin:$PATH
