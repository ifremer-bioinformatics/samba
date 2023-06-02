FROM amd64/ubuntu:20.04
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing R v4.2.3 and all R packages required for the samba pipeline"

RUN apt-get -y update
RUN apt-get -y install wget
RUN apt-get -y install build-essential
RUN apt-get -y install locales

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="$PATH:/opt/conda/bin"

RUN conda create -n R-4.2.3 -c conda-forge r-base procps-ng

SHELL ["conda", "run", "-n", "R-4.2.3", "/bin/bash", "-c"]
RUN conda install -c conda-forge mamba
RUN mamba install -c conda-forge glib
RUN mamba install -c conda-forge r-base
RUN mamba install -c conda-forge r-jpeg
RUN mamba install -c conda-forge r-futile.logger
RUN mamba install -c bioconda bioconductor-deseq2
RUN mamba install -c conda-forge r-stringr
RUN mamba install -c conda-forge r-dplyr
RUN mamba install -c conda-forge r-vegan
RUN mamba install -c conda-forge r-biocmanager
RUN R --vanilla -e 'BiocManager::install("phyloseq")'
RUN mamba install -c conda-forge r-tidyr
RUN mamba install -c conda-forge r-ggplot2
RUN mamba install -c conda-forge r-ggpubr
RUN mamba install -c conda-forge r-rstatix
RUN R --vanilla -e 'install.packages("grid", repos="https://mirror.ibcp.fr/pub/CRAN/")'
RUN mamba install -c conda-forge r-svglite
RUN mamba install -c conda-forge r-remotes
RUN R --vanilla -e 'remotes::install_github("mikemc/speedyseq")'

RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

RUN mamba install -c conda-forge r-reshape2
RUN mamba install -c r r-rcolorbrewer
RUN mamba install -c conda-forge r-lemon
RUN mamba install -c r r-rmarkdown
RUN mamba install -c conda-forge pandoc
RUN mamba install -c bioconda bioconductor-metagenomeseq
RUN mamba install -c conda-forge r-upsetr
RUN mamba install -c conda-forge r-plotly
RUN mamba install -c r r-htmlwidgets
RUN R --vanilla -e 'BiocManager::install("microbiome")'
RUN R --vanilla -e 'remotes::install_github("vmikk/metagMisc")'


RUN conda clean -a
ENV PATH /opt/conda/envs/R-4.2.3/bin:$PATH
