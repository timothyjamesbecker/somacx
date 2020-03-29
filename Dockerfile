FROM ubuntu:bionic
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
build-essential \
g++ \
gfortran \
git \
wget \
nano \
libffi6 \
libffi-dev \
libssl1.0.0 \
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
RUN python3 -m pip install cython
RUN python3 -m pip install numpy
RUN python3 -m pip install pysam
RUN python3 -m pip install https://github.com/timothyjamesbecker/somacx/releases/download/0.1.2/somacx-0.1.2.tar.gz
