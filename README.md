# DUESSELPORE Webserver manual

### Getting start
Simply install Docker and run:
```console
sudo docker run -it -p 8000:8000 thachdt4/duesselpore:running python3 /home/ag-rossi/projects/duesselpore/manage.py runserver 0.0.0.0:8000
```
Then open the web browser at localhost:8000/duesselpore. More information is on manual.md file.
Duesselpore is a full stack, light weigh webserver for Nanopore RNAseq.
<<<<<<< HEAD
* Duesselpore web framework is writen using Django Python 3.8.5, NGS data analysis function were written using R and Python.
* Duesselpore can run locally on your Docker Container.
* Most of functional analysis can work without internet connection.
=======
* Duesselpore web framework is writen using Django Python 3.8.5, NGS data analysis function were written using R and Python
* Duesselpore can run locally on your Docker Container or virtual machine (Virtualbox)
* After installation step (docker pull) most of functional analysis can work without internet connection.<br>

If you have any question, please contact: <br>
Thach Nguyen<br>
Genome Engineering and Model Development lab (GEMD)<br>
thach.nguyen(remove this bracket &%&//) at IUF-Duesseldorf.de

>>>>>>> 7070e16a91280dded774ddc30b9affddec6e39d3

### Availability:
* Running Docker images is on dockerhub: thachdt4/duesselpore:main, or you can customize your Docker container from Dockerfile on this github repository.
* Test data full at https://iufduesseldorf-my.sharepoint.com/:u:/g/personal/thach_nguyen_iuf-duesseldorf_de/EWIk4CLauThHk61_5rItjEcBOP4CJstbyCN9yN3ty36A7g?e=zRUf1T
* Test data lightweight at https://iufduesseldorf-my.sharepoint.com/:u:/g/personal/thach_nguyen_iuf-duesseldorf_de/ES4BsdfJSKNHl-mDUR3BogcBEmdOawVTRy-eRXU3-XeG2A?e=Kq9O2e 
<<<<<<< HEAD
=======
* Virtualbox image at https://iufduesseldorf-my.sharepoint.com/:u:/g/personal/thach_nguyen_iuf-duesseldorf_de/EZB0I5s_Gq5MnPbg1g69WvsBACVULQFQ3s2wJjc-pyN3PA?e=Jh7Hwv 
* Instruction video (for Virtualbox configuration, not recommended) https://www.youtube.com/watch?v=-MvffifrQ1Q 
>>>>>>> 7070e16a91280dded774ddc30b9affddec6e39d3

### System requirement
#### System requirement
    * CPU: 2.0 GHz (64bits) 2 cores or higher
    * Memory: 8 GB or higher
    * Diskdrive: 100 GB (or higher) free space
    * Window 10, Linux (64 bits) or Mac

### Installation and manual
See manual.pdf (manual.md) file or our instruction video

### Dependancies

Duesselpore use all solftware built in docker images.
* R version 4.1.1, Bioconductor and its packages
Full R packages and session info in R.txt
* Python 3.8.5, required packages in requirements.txt
* minimap2 version 2.20-r1061 and samtools version 1.7
* HTSeq 0.13.5
* Salmon v1.5.2
<<<<<<< HEAD
=======

### Reference 
Duesselpore: a full-stack local web server for rapid and simple analysis of Oxford Nanopore Sequencing data
Christian Vogeley, Thach Nguyen, Selina Woeste, Jean Krutmann, Thomas Haarmann-Stemmann, Andrea Rossi
bioRxiv 2021.11.15.468670; doi: https://doi.org/10.1101/2021.11.15.468670
>>>>>>> 7070e16a91280dded774ddc30b9affddec6e39d3
