# Statement of tutorial objectives

The aim of this tutorial is to demonstrate a workflow for long-read differential gene expression analysis based on cDNA sequence data. This workflow is suitable for fastq sequence collections with a paired design (e.g. tumour/normal) where a reference genome sequence and its gene annotation are available. 

The tutorial is packaged with example data, so that the workflow can be replicated to address questions such as

* which genes are expressed in my study?
* which genes are upregulated in this tumour sample?
* show expression levels for gene *ENSG00000117523*

Editing of the workflow's configuration file, **`config.yaml`**, will allow the workflow to be run with different starting cDNA sequence collections, reference genomes, and with different statistical requirements for the selection of candiate genes.

## Methods utilised include: 

* **`conda`** for management of bioinformatics software installations
* **`snakemake`** for managing the bioinformatics workflow
* **`minimap2`** for mapping sequence reads to reference genome
* **`samtools`** for SAM/BAM handling and mapping statistics
* **`RSubread`** and **`DESeq2`** for differential expression analysis

## The computational requirements include: 

* Computer running Linux (Centos7, Ubuntu 18_10, Fedora 29) - also validated on macOS (Mojave) 
* At least 16 Gb RAM - swap space of at least 8Gb ideal
* At least 15 Gb spare disk space for analysis and indices
* Runtime with provided example data - approximately 20 minutes

\pagebreak

# Software installation

1. Most software dependencies are managed though **`conda`**. Install as described at  <br> [https://conda.io/docs/install/quick.html](https://conda.io/docs/install/quick.html). You will need to accept the license agreement during installation and we recommend that you allow the Conda installer to prepend its path to your `.bashrc` file when asked.
```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    bash
```
2. Download Nanopore tutorial & example files into a folder named `cDNA_DESeq2`. This tutorial requires the **`git-lfs`** large file support capabilities, which should be installed through **`Conda`** first
```
    conda install -c conda-forge git-lfs
    git lfs install
    git clone https://github.com/nanoporetech/ont_tutorial_cdna_deseq2.git cDNA_DESeq2
```
3. Change your working directory into the new `cDNA_DESeq2` folder 
```
    cd cDNA_DESeq2
```
4. Install Conda software dependencies
```
    conda env create --name cDNA_DESeq2 --file environment.yaml
```
5. Initialise the Conda environment 
```
    source activate cDNA_DESeq2
```


\pagebreak

# Introduction

Differential gene expression analysis aims to identify genes that show statistically (and magnitudinally) altered expression patterns in a biological system of interest, when compared to a biological control. The results of the differential expression analysis are presented in a quantitative form and therefore the degree of change (up or down regulation) can be ascertained for each gene identified.  
  
Differential gene expression analysis requires a "snapshot" of gene activity. In this context, gene activity corresponds to the quantity of messenger RNAs (mRNA) transcribed from each gene within the organism/tissue/culture being investigated. The greater the number of mRNA molecules observed from a given gene, the higher the expression level. In order to determine expression levels across the whole genome, sequence data specifically targeting the mRNA molecules can be generated. 

[Oxford Nanopore Technologies](https://nanoporetech.com) provides a number of [sequencing solutions](https://nanoporetech.com/rna) to allow users to generate a snapshot of gene expression. This can be achieved by both sequencing the mRNA [directly](https://store.nanoporetech.com/catalog/product/view/id/167/s/direct-rna-sequencing-kit/category/28/), or via a complementary DNA ([cDNA](https://store.nanoporetech.com/catalog/product/view/id/177/s/cdna-pcr/category/28/)) proxy. In contrast to short read sequencing technologies, entire mRNA transcripts can be captured as single reads. The example data provided with this tutorial is from a study based on the [PCR-cDNA](https://store.nanoporetech.com/catalog/product/view/id/177/s/cdna-pcr/category/28/) kit. This is a robust choice for performing differential gene expression studies. This kit is suitable for preparation of sequence libraries from low mRNA input quantities. The cDNA population is [enriched through PCR with low bias](https://nanoporetech.com/resource-centre/low-bias-rna-seq-pcr-cdna-pcr-free-direct-cdna-and-direct-rna-sequencing); an important prerequisite for the subsequent statistical analysis.   
  
Once sequencing data has been produced from both the experimental and paired control samples (with an appropriate number of biological replicates), the sequence reads can be mapped to the organism's reference genome. The number of sequences mapping to each gene can be counted, and it is this count data that forms the basis for differential gene expression analysis.  

There are five goals for this tutorial:

* To introduce a literate framework for analysing Oxford Nanopore cDNA data prepared using the MinION, GridION or PromethION
* To utilise best data-management practices
* To provide basic cDNA sequence QC metrics, enabling review and consideration of the starting experimental data
* To map sequence reads to the reference genome and to identify the genes that are expressed and the number of sequence reads that are observed from each gene
* To perform a statistical analysis using **`DESeq2`** (@R-DESeq2) to identify differentially expressed genes

This tutorial does not aim to provide an exhaustive analysis or annotation of the differentially expressed genes. 

# Getting started and best practices

This tutorial requires a computer workstation running a Linux operating system. The workflow described has been tested using **`Fedora 29`**, **`Centos 7`** and **`Ubuntu 18_04`**. This tutorial has been prepared in the **`Rmarkdown`** file format. This utilises *markdown* (an easy-to-write plain text format as used in many Wiki systems) - see @R-rmarkdown for more information about **`rmarkdown`**. The document template contains chunks of embedded **`R code`** that are dynamically executed during the report preparation. 

The described analytical workflow makes extensive use of the **`conda`** package management and the **`snakemake`** workflow software. These software packages and the functionality of **`Rmarkdown`** provide the source for a rich, reproducible and extensible tutorial document.

The workflow contained within this Tutorial performs an authentic bioinformatics analysis and using the whole human genome as a reference sequence. There are some considerations in terms of memory and processor requirement. Indexing the whole human genome for sequence read mapping using **`minimap2`** for example will use at least **`18 Gb`** of memory. The minimal recommended hardware setup for this tutorial is therefore an 4 threaded computer with at least 16 Gb of RAM and 15 Gb of storage space. 

There are few dependencies that need to be installed at the system level prior to running the tutorial. The **`conda`** package management software will coordinate the installation of the required bioinformatics software and their dependencies in user space - this is dependent on a robust internet connection.

As a best practice this tutorial will separate primary DNA sequence data (the base-called fastq files) from the **`Rmarkdown`** source and the genome reference data. The analysis results and figures will again be placed in a separate working directory. The required layout for the primary data is shown in the figure below. This minimal structure will be prepared over the next few sections of this tutorial. The DNA sequences must be placed within a folder called **`RawData`** and the reference genome and annotation files must be placed in a folder named **`ReferenceData`**.


![](Static/Images/FolderLayout.png) 

# Experimental setup

The first required step for performing a cDNA differential sequence analysis involves collation of information on the biological samples and biological replicates that are to be used for the statistical analysis.

![](Static/Images/ExperimentalDesign.png) 

The example data included with this tutorial describes a study comparing an experimental sample against a linked control. The experimental and control samples have been prepared in triplicate. This design is described in a configuration file named **`config.yaml`** - an example file has been provided with the tutorial. The content of this file is highlighted in the figure above. The cDNA sequence files are defined within the **`Sample`** block; experimental groups and their discrete biological samples are defined within this block.

**`reference_genome`** refers to the genome against which the sequence reads will be mapped; **`genome_annotation`** refers to the gene annotations assigned to this genome sequence. In this tutorial a URL is provided for both and the **`Snakemake`** workflow will download the corresponding files. 

The only other parameters that should be considered are **`readCountMinThreshold`**, **`lfcThreshold`** and **`adjPValueThreshold`**. These correspond to the required depth of coverage to analyse a gene, the log2-fold-change filter to be applied in differential testing and the false-discovery corrected p-value threshold to be used, respectively.


## Example dataset

This tutorial is distributed with a collection of Oxford Nanopore cDNA sequence data in fastq format. A sequence collection corresponding to a renal cancer sample and its corresponding normal have been prepared. The sequence collection has been filtered and sub-sampling has been performed to yield a **synthetic dataset** that can be used to demonstrate this workflow. 

\newpage

# Snakemake

This tutorial for expression profiling from cDNA sequence data makes use of **`snakemake`** (@snakemake2012). Snakemake is a workflow management system implemented in Python. The aim of the Snakemake tool is to enable reproducible and scalable data analyses. The workflow produced within this document should be portable between laptop resources, computer servers and other larger scale IT deployments. The Snakemake workflow additionally defines the sets of required software (and software versions where appropriate) and will automate the installation and deployment of this software through the **conda** package management system.

The **`snakemake`** workflow will call methods that include **`minimap2`** @minimap22018 and **`samtools`** @samtools2009. The planned workflow is shown in the figure below. The provided reference genome sequence will be indexed using **`minimap2`**, each sequence collection will be mapped to the index (again using **`minimap2`** with parameters tuned for the mapping of long reads whilst accommodating exon matches interspersed by introns) and summary statistics will be prepared using the **`samtools`** software. The remainder of the analysis will be performed in the **`R analysis`** described within the report.

![](Static/Images/dag1.png) 

The precise commands within the **`Snakefile`** include

* download the specified reference genome
* download the specified genome annotations
* use **`minimap2`** to index the reference genome
* map cDNA sequence reads against the reference genome index using **`minimap2`**
* convert **`minimap2`** output (**`SAM`**) into a sorted **`BAM`** format using **`samtools`**
* prepare summary mapping statistics using **`samtools`**

# Run the snakemake workflow file

The snakemake command is responsible for orchestrating the analytical workflow. The snakemake command with the parameters **`snakemake --forceall --dag | dot -Tpng | display`** will prepare and display the a flowchart that describes the workflow. The flowchart for this tutorial is displayed in the figure above.

\fontsize{8}{12}
```
# just type snakemake to run the workflow
# don't type <NPROC> but specify the number of processor cores available (e.g. 2 or 4)

snakemake -j <NPROC>
```
\fontsize{10}{14}


\pagebreak
 

# Prepare the analysis report

The **`Rmarkdown`** script can be run usimg the **`knit`** dialog in the **`Rstudio`** software - please see the figure below. Selecting **`Knit to HTML`** will prepare a portable HTML file. 

![](Static/Images/KnitIt.png) 

The document can also be rendered from the command line with the following command

\fontsize{8}{12}
```
R --slave -e 'rmarkdown::render("Nanopore_cDNA_Tutorial.Rmd", "html_document")'
```
\fontsize{10}{14}

