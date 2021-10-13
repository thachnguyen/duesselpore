FROM bioconductor/bioconductor_docker:devel

# File Author / Maintainer
MAINTAINER Thach Nguyen <thach.nguyen@iuf-dueseldorf.de> 

ARG DUESSELPORE=1.0.0

# setup working directory
RUN mkdir /home/ag-rossi; cd /home/ag-rossi; mkdir projects

#Install Salmon 
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz; tar xzvf salmon-1.5.2_linux_x86_64.tar.gz; rm salmon-1.5.2_linux_x86_64.tar.gz 
RUN ln -s /home/ag-rossi/salmon-1.5.2_linux_x86_64/bin/salmon /bin/salmon

#Install minimap2
RUN git clone https://github.com/lh3/minimap2; cd minimap2 && make; ln -s /home/ag-rossi/minimap2/minimap2 /bin/minimap2

RUN cd /home/ag-rossi/projects/; git clone https://github.com/thachnguyen/duesselpore.git; chmod +x duesselpore 
#Install R packages
RUN cd duesselpore; R < R.txt --no-save

# Install python dependancies
RUN pip install django==3.2.8 requests==2.26.0 biopython==1.78 HTSeq==0.13.5 openpyxl==3.0.9 pyensembl==1.9.4

#Install decompressor
RUN apt-get update; apt-get -y install p7zip-full; apt install -y fastqc; apt-get install -y samtools

#Fixed fastqc bug 
RUN cp /home/ag-rossi/projects/duesselpore/duesselpore/fastqc /usr/bin/fastqc
RUN python3 

RUN cd /home/ag-rossi/projects/duesselpore; python3 setup_webserver_docker.py light
# if you want to use other orgarnisms, run this command  
# python3 setup_webserver_docker.py full 
RUN python3 manage.py runserver 0.0.0.0:8000
ENV LC_ALL=en_US.utf8
ENV LANG=en_US.utf8
