FROM amd64/ubuntu:20.04
LABEL authors="samba-sebimer@ifremer.fr" \
      description="Docker image containing FIGARO-1.1.2 requirements for the samba pipeline"

RUN apt -y update
RUN apt -y install git
RUN apt -y install python3
RUN apt -y install python3-pip
RUN pip3 install setuptools
RUN mkdir -p /opt
WORKDIR /opt
RUN git clone https://github.com/Zymo-Research/figaro.git
WORKDIR /opt/figaro
RUN python3 setup.py bdist_wheel
RUN pip3 install --force-reinstall dist/*.whl
