FROM ubuntu:focal-20250404
LABEL org.opencontainers.image.authors="timothyjamesbecker@gmail.com"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
build-essential \
g++ \
gfortran \
git \
wget \
nano \
libffi7 \
libffi-dev \
libssl1.1 \
libssl-dev \
libblas3 \
libblas-dev \
liblapack3 \
liblapack-dev \
libcurl4-openssl-dev \
libxml2-dev \
libncurses5-dev \
libbz2-dev \
zlib1g-dev \
python3 \
python3-dev \
python3-pip
RUN apt-get autoremove
RUN python3 -m pip install pip --upgrade pip
RUN python3 -m pip install -Iv cython==3.0.12 #latest working build versions here
RUN python3 -m pip install -Iv numpy==1.24.4  #latest working build versions here
RUN python3 -m pip install -Iv pysam==0.23.0  #latest working build versions here
RUN python3 -m pip install https://github.com/timothyjamesbecker/somacx/releases/download/0.1.3/somacx-0.1.3.tar.gz
RUN unbuffered generator.py -h && echo "somacx-0.1.3 install was successful!"
