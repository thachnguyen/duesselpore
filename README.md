# DUESSELPORE Webserver manual

### Getting start
Simply install Docker and run:
```console
sudo docker run -it -p 8000:8000 thachdt4/duesselpore:running python3 /home/ag-rossi/projects/duesselpore/manage.py runserver 0.0.0.0:8000
```
Then open the web browser at localhost:8000/duesselpore. More information is on manual.md file.
Duesselpore is a full stack, light weigh webserver for Nanopore RNAseq.
* Duesselpore web framework is writen using Django Python 3.8.5, NGS data analysis function were written using R and Python
* Duesselpore can run locally on your Docker Container or virtual machine (Virtualbox)
* After installation step (docker pull) most of functional analysis can work without internet connection.

If you have any question, please contact: 
Thach Nguyen
Genome Engineering and Model Development lab (GEMD)
thach.nguyen(remove this bracket &%&//) at IUF-Duesseldorf.de


### Availability:
* Running Docker images is on dockerhub: thachdt4/duesselpore:main, or you can build your Docker container from Dockerfile on this github repository.    
* Test data full at https://iufduesseldorf-my.sharepoint.com/:u:/g/personal/thach_nguyen_iuf-duesseldorf_de/EWIk4CLauThHk61_5rItjEcBOP4CJstbyCN9yN3ty36A7g?e=zRUf1T
* Test data lightweight at https://iufduesseldorf-my.sharepoint.com/:u:/g/personal/thach_nguyen_iuf-duesseldorf_de/ES4BsdfJSKNHl-mDUR3BogcBEmdOawVTRy-eRXU3-XeG2A?e=Kq9O2e 
* Virtualbox image at https://iufduesseldorf-my.sharepoint.com/:u:/g/personal/thach_nguyen_iuf-duesseldorf_de/EZB0I5s_Gq5MnPbg1g69WvsBACVULQFQ3s2wJjc-pyN3PA?e=Jh7Hwv 
* Instruction video (for Virtualbox configuration, not recommended) https://www.youtube.com/watch?v=-MvffifrQ1Q 

### System requirement
#### Minimum requirement
    * CPU: 2.0 GHz (64bits CPU acceleration support) 4 cores or higher
    * Memory: 8 GB or higher
    * Diskdrive: 200 GB free space
    * Window 10, Linux (Ubuntu) or Mac

#### Recommend configuration

    * CPU: 3.0 GHz (64bits) 8 cores or higher
    * Memory: 16 GB or higher
    * Diskdrive: 1000 GB free space
    * Window 10, Linux (Ubuntu) or Mac (64 Operating system)

### Installation and manual
See manual.pdf (manual.md) file or our instruction video

### Dependancies

Duesselpore use 
* R version 4.1.1, Bioconductor and its packages
Full R packages and session info in R.txt
* Python 3.8.5, required packages in requirements.txt
* minimap2 version 2.20-r1061 and samtools version 1.7
* HTSeq 0.13.5
* Salmon v1.5.2

### Reference 
Duesselpore: a full-stack local web server for rapid and simple analysis of Oxford Nanopore Sequencing data
Christian Vogeley, Thach Nguyen, Selina Woeste, Jean Krutmann, Thomas Haarmann-Stemmann, Andrea Rossi
bioRxiv 2021.11.15.468670; doi: https://doi.org/10.1101/2021.11.15.468670
