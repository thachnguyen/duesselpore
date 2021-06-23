# NGS_webserver the webserver for RNA Seq data analysis

# Fastq file is ONE zip file contains all the replicates of each experimental group. The directory name is also the group name.  
Examples the directory structure showed below, the replicate's file name are arbitrary:
mydata.zip/treated01/replica01.fastq
	            /replica02.fastq
	  /untreated01/replica03.fastq
		      /replica04.fastq
If one of your group contain only one replica. The statistical approach will be disable and the interception approach base on pairwise comparison will activate.

# Gene count matrices (optional): Instead of processing fastq files: our webserver support precomputed gene count matrices

# Reference genome: Our webserver support both official reference genome and customize reference genome as fasta format 

# The parameter of DESeq2 parameter:
- ReadCountMinThreshold (Optional): the Minimum read count
- LfcThreshold (Optional): log2 threshold
- AdjPValueThreshold (Optional): adjudgement P value threshold 




