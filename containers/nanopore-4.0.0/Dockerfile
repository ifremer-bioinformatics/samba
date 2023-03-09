FROM --platform=linux/x86_64 mambaorg/micromamba
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing all dependancies required for the Nanopore samba pipeline"

COPY environment.yml /
RUN micromamba create -f /environment.yml
RUN micromamba clean -a
ENV PATH /opt/conda/envs/nanopore-4.0.0/bin:$PATH
