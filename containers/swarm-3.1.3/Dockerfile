FROM --platform=linux/x86_64 mambaorg/micromamba
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing swarm 3.1.3 required for the samba pipeline"

COPY environment.yml /
RUN micromamba create -f /environment.yml
RUN micromamba clean -a
ENV PATH /opt/conda/envs/swarm-3.1.3/bin:$PATH
