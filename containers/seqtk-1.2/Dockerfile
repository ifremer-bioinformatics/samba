FROM nfcore/base:1.9
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing seqtk requirements for the nfcore/samba pipeline"

COPY containers/seqtk-1.2/environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-seqtk-1.2/bin:$PATH
