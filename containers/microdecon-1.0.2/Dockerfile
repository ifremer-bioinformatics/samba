FROM nfcore/base:1.9
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing r-microDecon requirements for the nfcore/samba pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-microdecon-1.0.2/bin:$PATH
