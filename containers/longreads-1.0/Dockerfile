FROM nfcore/base:1.9
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing long reads requirements for the samba pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/longreads-1.0/bin:$PATH
