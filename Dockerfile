# Base Image
FROM immcantation/suite:4.4.0

LABEL maintainer="VDJServer <vdjserver@utsouthwestern.edu>"

# extra tools
RUN yum install -y jq

# AIRR v1.6.0
RUN pip3 uninstall -y airr
#RUN pip3 install airr
RUN git clone https://github.com/airr-community/airr-standards.git
RUN cd airr-standards && git checkout v1.6.0 && cd lang/python && pip3 install .

# Copy source
RUN mkdir /repcalc-root
COPY . /repcalc-root

# install repcalc
RUN cd /repcalc-root && pip3 install .
