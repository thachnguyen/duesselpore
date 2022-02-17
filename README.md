# DUESSELPORE Webserver manual

### Getting start
Simply install Docker and run (you may have to run it as sudo user on Linux):
```console
docker run -it -p 8000:8000 thachdt4/duesselpore:running python3 /home/ag-rossi/projects/duesselpore/manage.py runserver 0.0.0.0:8000
```
Open the web browser at localhost:8000/duesselpore. More information is on manual.md file.
Duesselpore is a full stack, light weight webserver to process RNASeq from Nanopore as well as other (Illumina, PacBio ...)
* Duesselpore web framework is writen using Django Python 3.8.5, NGS data analysis function were written using R and Python.
* Duesselpore can run locally on your Docker Container.
* Most of functional analysis can work without internet connection.

### Availability:
* Running Docker images is on dockerhub: thachdt4/duesselpore:main, or you can customize your Docker container from Dockerfile on this github repository.
* Test data full at https://iufduesseldorf-my.sharepoint.com/:u:/g/personal/thach_nguyen_iuf-duesseldorf_de/EWIk4CLauThHk61_5rItjEcBOP4CJstbyCN9yN3ty36A7g?e=zRUf1T (Human)
* Test data lightweight at https://iufduesseldorf-my.sharepoint.com/:u:/g/personal/thach_nguyen_iuf-duesseldorf_de/ES4BsdfJSKNHl-mDUR3BogcBEmdOawVTRy-eRXU3-XeG2A?e=Kq9O2e (Human)
* Instruction video 

### System requirement
#### System requirement
    * CPU: 2.0 GHz (64bits) 2 cores or higher
    * Memory: 8 GB or higher
    * Diskdrive: 100 GB free space 
    * Window 10, Linux (64 bits) or MacOS.

### Dependancies
All dependancies package are built in docker images.
