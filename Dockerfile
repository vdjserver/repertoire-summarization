# Base Image
FROM immcantation/suite:4.4.0

LABEL maintainer="VDJServer <vdjserver@utsouthwestern.edu>"

# extra tools
RUN yum install -y jq

# latest AIRR
RUN pip3 uninstall -y airr
RUN pip3 install airr

# Copy source
RUN mkdir /repcalc-root
COPY . /repcalc-root

# install repcalc
RUN cd /repcalc-root && pip3 install .
