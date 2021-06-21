# Base Image
FROM ubuntu:20.04

MAINTAINER VDJServer <vdjserver@utsouthwestern.edu>

# uncomment these if building behind UTSW proxy
#ENV http_proxy 'http://proxy.swmed.edu:3128/'
#ENV https_proxy 'https://proxy.swmed.edu:3128/'
#ENV HTTP_PROXY 'http://proxy.swmed.edu:3128/'
#ENV HTTPS_PROXY 'https://proxy.swmed.edu:3128/'

# Install OS Dependencies
RUN export DEBIAN_FRONTEND=noninteractive && apt-get update && apt-get install -y \
    libssl-dev \
    python3 \
    python3-pip \
    python3-scipy \
    r-base \
    r-base-dev \
    libssh2-1-dev \
    libcurl4-openssl-dev \
    libyaml-dev \
    mercurial

RUN pip install \
    numpy \
    lxml \
    argparse \
    BeautifulSoup4 \
    reportlab \
    biopython \
    airr

RUN pip3 install \
    pandas \
    biopython \
    airr \
    presto \
    changeo \
    sphinx \
    sphinxcontrib-autoprogram \
    prov \
    recommonmark \
    rdflib \
    networkx \
    lxml \
    commonmark \
    isodate \
    decorator \
    future \
    agavepy

# Copy source
RUN mkdir /repcalc-root
COPY . /repcalc-root

# install repcalc
RUN cd /repcalc-root && python3 setup.py install
