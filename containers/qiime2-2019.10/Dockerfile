FROM nfcore/base:1.9
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing QIIME2-2019.10 requirements for the nfcore/samba pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-qiime2-2019.10/bin:$PATH
