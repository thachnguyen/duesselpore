FROM bioconductor/bioconductor_docker:devel

# File Author / Maintainer
MAINTAINER Thach Nguyen <thach.nguyen@iuf-dueseldorf.de> 

ARG DUESSELPORE=1.0.0

#Install minimap2, setup working directory
RUN mkdir /home/ag-rossi; cd /home/ag-rossi; mkdir projects
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz; tar xzvf salmon-1.5.2_linux_x86_64.tar.gz; rm salmon-1.5.2_linux_x86_64.tar.gz 

RUN git clone https://github.com/lh3/minimap2; cd minimap2 && make; ln -s /home/ag-rossi/minimap2/minimap2 /bin/minimap2

RUN cd /home/ag-rossi/projects/; git clone https://github.com/thachnguyen/duesselpore.git 
#Install R packages
RUN cd /home/ag-rossi/projects/duesselpore
RUN R < R.txt --no-save

#Install Salmon 
RUN ln -s /home/ag-rossi/salmon-1.5.2_linux_x86_64/bin/salmon /bin/salmon

# Install python dependancies
RUN pip install django==3.2.8 requests==2.26.0 biopython==1.78 HTSeq==0.13.5

#Install decompressor
RUN apt-get update
RUN apt-get install -y p7zip-full

RUN cd /home/ag-rossi/projects/duesselpore; python3 manage.py runserver 0.0.0.0:8000
ENV LC_ALL=en_US.utf8
ENV LANG=en_US.utf8
