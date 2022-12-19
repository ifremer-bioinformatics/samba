FROM --platform=linux/x86_64 mambaorg/micromamba
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing QIIME2-2022.8 requirements for the samba pipeline"

COPY environment.yml /
RUN micromamba create -f /environment.yml && micromamba clean -a
ENV PATH /opt/conda/envs/qiime2-2022.8/bin:$PATH
