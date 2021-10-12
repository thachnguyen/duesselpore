FROM bioconductor/bioconductor_docker:devel

# File Author / Maintainer
MAINTAINER Thach Nguyen <thach.nguyen@iuf-dueseldorf.de> 

ARG DUESSELPORE=1.0.0

RUN mkdir /home/ag-rossi; cd ag-rossi; mkdir projects
RUN git clone https://github.com/lh3/minimap2; cd minimap2 && make; ln -s /home/ag-rossi/minimap2/minimap2 /bin/minimap2

RUN git clone https://github.com/thachnguyen/duesselpore.git 
RUN python3 

RUN pip install django==3.2.8 requests==2.26.0 biopython==1.78
RUN apt-get update; apt-get install p7zip-full

# Installing modules for EPINANO
RUN pip install h5py==2.8.0 numpy==1.16.1 pandas==0.23.4 scikit-learn==0.20.2
RUN yum install -y wget 
RUN wget https://github.com/enovoa/EpiNano/archive/refs/tags/epinano${EPINANO_VERSION}.tar.gz
RUN tar -zvxf epinano${EPINANO_VERSION}.tar.gz; cd EpiNano-epinano${EPINANO_VERSION}/
RUN chmod +x /project/EpiNano-epinano${EPINANO_VERSION}/scripts/*
ENV PATH "/project/EpiNano-epinano${EPINANO_VERSION}/scripts/" 

# Clean cache
RUN yum clean all 

#cleaning
RUN rm -fr *.tar.gz; rm -fr *.bz2; rm -fr ./tmp
RUN rm -rf /var/cache/yum

ENV LC_ALL=en_US.utf8
ENV LANG=en_US.utf8
