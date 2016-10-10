# Base Image
FROM debian:jessie

MAINTAINER VDJServer <vdjserver@utsouthwestern.edu>

# uncomment these if building behind UTSW proxy
ENV http_proxy 'http://proxy.swmed.edu:3128/'
ENV https_proxy 'https://proxy.swmed.edu:3128/'
ENV HTTP_PROXY 'http://proxy.swmed.edu:3128/'
ENV HTTPS_PROXY 'https://proxy.swmed.edu:3128/'

# Install OS Dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    doxygen \
    git \
    graphviz \
    libbz2-dev \
    libxml2-dev \
    libxslt-dev \
    python \
    python-dev \
    python-sphinx \
    python-pip \
    vim \
    wget \
    zlib1g-dev \
    cpio \
    emacs

RUN pip install \
    biopython \
    lxml \
    numpy \
    argparse \
    BeautifulSoup4 \
    reportlab

# Boost
ENV BOOST_VERSION 1.57.0
ENV BOOST_VERSION_LINK 1_57_0

# Install/bootstrap boost
RUN wget http://downloads.sourceforge.net/project/boost/boost/$BOOST_VERSION/boost_$BOOST_VERSION_LINK.tar.gz
RUN tar -xvzf boost_$BOOST_VERSION_LINK.tar.gz
RUN cd /boost_$BOOST_VERSION_LINK && ./bootstrap.sh --prefix=/usr/local
RUN cd /boost_$BOOST_VERSION_LINK && ./b2 install
RUN cd /boost_$BOOST_VERSION_LINK/tools/build && ./bootstrap.sh
RUN cd /boost_$BOOST_VERSION_LINK/tools/build && ./b2 install --prefix=/usr/local

# VDJML
ENV VDJML_VERSION 0.1.4
RUN git clone https://bitbucket.org/vdjserver/vdjml.git vdjml-root
RUN cd /vdjml-root && git checkout develop
RUN cp /vdjml-root/docker/boost/boost-build.jam /
RUN cp /vdjml-root/docker/boost/user-config.jam /root/

RUN cd /vdjml-root && b2
RUN cd /vdjml-root && b2 distro-bindings-py
RUN tar zxvf /vdjml-root/out/VDJMLpy-$VDJML_VERSION.tar.gz

# extract database
COPY db.tgz /
RUN tar zxvf db.tgz

# Copy source
RUN mkdir /repsum-root
COPY . /repsum-root

# install repsum
RUN cd /repsum-root && python setup.py install

# setup run environment
ENV DB_DIR "/db"
ENV VDJ_DB_ROOT "/db/10_05_2016/"
ENV IGDATA "/db"
ENV VDJSERVER_ROOT "/repsum-root"
ENV PATH "$PATH:/igblast-root/local/bin"
ENV PYTHONPATH "/VDJMLpy-$VDJML_VERSION:$PYTHONPATH"
ENV PYTHONPATH "/VDJMLpy-$VDJML_VERSION/vdjml:$PYTHONPATH"
ENV LD_LIBRARY_PATH "/VDJMLpy-$VDJML_VERSION/vdjml:$LD_LIBRARY_PATH"
