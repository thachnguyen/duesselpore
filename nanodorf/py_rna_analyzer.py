'''
This is the main application functions. Read the user's files (fastaq format). Align to the reference genome use Minimap2. Analyse the map to the reference genome.
Export the read and Gene mapping result to an excel file. 
'''
from django.apps import AppConfig
import glob
from Bio import SeqIO
import os
from zipfile import ZipFile
import shutil
import zipfile
from django.core.mail import send_mail

from django.http import HttpResponse
import numpy as np


def read_fastq_list(s_id):
    '''
    Read user upload files in zip format and return the list of fastq entries for later process
    '''
    #Note: This method ignore fastq files' name.
    zipfilename = 'users_file/%s/fastq.zip' %s_id 
    file1 = ZipFile(zipfilename)
    path = 'tmp/fastq/%s/'%s_id
    os.mkdir(path)
    file1.extractall(path=path)
    fastq_files = glob.glob(path+'*.fastq')
    sq_list = []
    for f in fastq_files:
        seq = list(SeqIO.parse(f, format='fastq'))
        sq_list += seq
    return sq_list   

# seq_list = read_fastq_list(s)

def run_minimap2(path='users_file/', s_id = 'Test_name_1618217069'):
    os.mkdir('users_file/%s/Analysis'%s_id)
    path_minimap = 'users_file/%s/Analysis/Minimap'%s_id
    os.mkdir(path_minimap)
    path_flagstat = 'users_file/%s/Analysis/flagstat/'%s_id
    os.mkdir(path_flagstat)
    path1 = '%s%s'%(path,s_id)
    for group in os.listdir('users_file/%s/fastq/'%s_id):
        for fastq_file in os.listdir('users_file/%s/fastq/%s'%(s_id, group)):
            path2 = '%s/fastq/%s/%s'%(path1, group, fastq_file)
            fastq_file1 = fastq_file.split('.')[0]
            print(('minimap2 -t16 -ax splice -k14 -secondary=no /home/ag-rossi/ReferenceData/reference.mmi %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(path2, path_minimap, fastq_file1)))
            print(('samtools flagstat %s/%s.bam>%s%s.txt'%(path_minimap, fastq_file1, path_flagstat, fastq_file1)))
            os.system('minimap2 -t16 -ax splice -k14 -secondary=no /home/ag-rossi/ReferenceData/reference.mmi %s | samtools view -Sb | samtools sort - -o %s/%s.bam'%(path2, path_minimap, fastq_file1))
            os.system('samtools flagstat %s/%s.bam>%s%s.txt'%(path_minimap, fastq_file1, path_flagstat, fastq_file1))
    return


def create_yaml(s_id, samples, yaml_file = 'config.yaml', ref_group = 0, readCountMinThreshold = 10, lfcThreshold =1,  adjPValueThreshold = 0.05, tutorialText=False):
    # '''
    # Create the user defined YAML from YAML template for RScript. The fastq files are separate into different groups as subfolders. Named by directory's name and file'sname.
    # '''
    import yaml
    with open(yaml_file) as file:
        default_config = yaml.load(file, Loader=yaml.FullLoader)
    default_config['readCountMinThreshold'] = readCountMinThreshold
    default_config['lfcThreshold'] = lfcThreshold
    default_config['adjPValueThreshold'] = adjPValueThreshold
    default_config['tutorialText'] = tutorialText
    default_config['Samples'] = samples
    with open('users_file/%s/config.yaml'%s_id, 'w') as f:
        yaml.dump(default_config, f)

    return

# def write_rscript(path=):
#     template_R =
#     'setwd(\\"/home/ag-rossi/projects/NGS_webserver/NGS_webserver/RNAseq_analyzer/%s/\\")\n\
#     library(digest)\n\
#     library(ShortRead)\n\
#     library(ggplot2)\n\
#     library(plyr)\n\
#     library(dplyr)\n\
#     library(tidyr)\n\
#     library(Rsubread)\n\
#     library(DESeq2)\n\
#     library(pcaMethods)\n\
#     library(kableExtra)\n\
#     library(caTools)\n\
#     library(writexl)\n\
#     library(yaml)\n\
#     library(session)\n\
#     library(AnnotationDbi)\n\
#     library(org.Hs.eg.db)\n\
#     source("Static/R/common.R")\n\
#     resultDir <- file.path("Analysis", "Results")\n\
#     dir.create(resultDir, showWarnings = FALSE, recursive = TRUE)\n\
#     config <- yaml::yaml.load_file("config.yaml")\n\
#     persistenceData <- file.path(resultDir, "NanoporeDESeq2.Rdata")\n
#     studyDesign <- data.frame() \nfor (i in 1:length(config$Samples)) {\n  studyDesign <- rbind(studyDesign,\n                       data.frame(samples=names(config$Samples[[i]][[1]]), \n                                  filename=unlist(config$Samples[[i]][[1]]), \n                                  group=names(config$Samples[[i]])))\n}\nstudyDesign$replicate <- sapply(1:nrow(studyDesign), function(x)sum(studyDesign$group[1:x]==studyDesign$group[x]))\nstudyDesign$md5 <- lapply(as.character(studyDesign$filename), md5sum)\nstudyDesign <- studyDesign[,-which(colnames(studyDesign)=="samples")]\nprocessQCFastq <- function(rowname) {\n  row <- which(row.names(studyDesign)==rowname)\n  file <- as.character(studyDesign[row, "filename"])\n  fastq <- readFastq(file, widthIds=FALSE)\n  c(\n    reads = formatC(length(fastq), big.mark=","),\n    mbs = formatC(round(sum(width(fastq)) / 1000 / 1000, digits=1), big.mark=","),\n    min = min(width(fastq)),\n    max = max(width(fastq)),\n    mean = round(mean(width(fastq)), digits=1),\n    median = round(median(width(fastq)), digits=0),\n    qval = round(mean(alphabetScore(fastq) / width(fastq)), digits=1),\n    gc = round(mean(letterFrequency(sread(fastq), "GC")  / width(fastq)) * 100, digits=1),\n    n50 = ncalc(width(fastq), n=0.5), # The ncalc function is defined in the common.R-File\n    l50 = lcalc(width(fastq), n=0.5),\n    n90 = ncalc(width(fastq), n=0.9),\n    l90 = lcalc(width(fastq), n=0.9)\n  )\n}\ndata <- lapply(row.names(studyDesign), processQCFastq)\nqcData <- data.frame(data)\ncolnames(qcData) <- row.names(studyDesign)\nView(qcData)\nknitr::kable(qcData, caption="Summary statistics for the cDNA libraries imported", booktabs=TRUE, table.envir=\'table*\', linesep="")  %>%\n  kable_styling(latex_options=c("hold_position", font_size=9))\nextractLengths <- function(rowname) {\n  row <- which(row.names(studyDesign)==rowname)\n  file <- as.character(studyDesign[row, "filename"])\n  fastq <- readFastq(file)\n  t(as.matrix(width(fastq)))\n}\nlengthData <- mapply(extractLengths, row.names(studyDesign))\nlengthDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(lengthData)))\ncolnames(lengthDataMatrix) <-  row.names(studyDesign)\n\nlengthMatrixMelt <- reshape2::melt(lengthDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))\nlengthMatrixMelt <- cbind(lengthMatrixMelt, group=studyDesign[match(lengthMatrixMelt$variable,  rownames(studyDesign)), "group"])\n\nplot <- ggplot(lengthMatrixMelt, aes(x=variable, y=value, fill=group)) + geom_violin() + scale_y_continuous(limits=c(0, as.numeric(quantile(lengthMatrixMelt$value, probs=c(0.975))))) + xlab("study sample") +  ylab("Distribution of Read Lengths (bp)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution of read lengths across samples")\n\nsuppressWarnings(print(plot))\n\nextractQualities <- function(rowname) {\n  row <- which(row.names(studyDesign)==rowname)\n  file <- as.character(studyDesign[row, "filename"])\n  fastq <- readFastq(file)\n  t(as.matrix(alphabetScore(fastq) / width(fastq)))\n}\nqualityData <- mapply(extractQualities, row.names(studyDesign))\nqualityDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(qualityData)))\ncolnames(qualityDataMatrix) <-  row.names(studyDesign)\n\nqualityMatrixMelt <- reshape2::melt(qualityDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))\nqualityMatrixMelt <- cbind(qualityMatrixMelt, group=studyDesign[match(qualityMatrixMelt$variable,  rownames(studyDesign)), "group"])\n\nplotQ <- ggplot(qualityMatrixMelt, aes(x=variable, y=value, fill=group)) + geom_violin() + scale_y_continuous(limits=c(min(qualityMatrixMelt$value), max(qualityMatrixMelt$value))) + xlab("study sample") +  ylab("Distribution of Read Qualities (QV)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution of read qualities across samples")\n\nsuppressWarnings(print(plotQ))\n\nflagstatTargets <- file.path("Analysis", "flagstat", \n                             paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE),".txt",sep=""))\n\nloadFlagstat <- function(file) {\n  x <- read.table(file, header=FALSE, sep=" ", fill=NA)[c(1:5),c(1,3)]\n  x[,1]\n}\n\nflagstatRes <- data.frame(matrix(unlist(lapply(flagstatTargets, loadFlagstat)), ncol=length(flagstatTargets)), stringsAsFactors = FALSE)\ncolnames(flagstatRes) <- rownames(studyDesign)\nrownames(flagstatRes) <- c("read mappings", "Secondary", "Supplementary", "Duplicates", "Mapped")\n\nflagstatRes[nrow(flagstatRes)+1,] <- as.numeric(gsub(",","",t(qcData)[, "reads"]))\nrownames(flagstatRes)[6] <- "nreads"\n\ngetVal <- function(word) {\n  sum(as.numeric(unlist(strsplit(word, "/"))))\n}\n\nzreads <- unlist(lapply(flagstatRes["read mappings", ], getVal)) -\n  unlist(lapply(flagstatRes["Secondary", ], getVal)) -\n  unlist(lapply(flagstatRes["Supplementary", ], getVal)) - \n  unlist(lapply(flagstatRes["Duplicates", ], getVal)) \n\nflagstatRes[nrow(flagstatRes)+1,] <- round(as.numeric(flagstatRes["Mapped", ]) / as.numeric(flagstatRes["nreads", ]) * 100, digits = 2)\n\nflagstatRes <- flagstatRes[c(6,1,2,3,4,7),]\nrownames(flagstatRes)[6] <- "%mapping"\n\nknitr::kable(flagstatRes, caption="Summary statistics from the minimap2 long read spliced mapping.", booktabs=TRUE, table.envir=\'table*\', linesep="")  %>%\n  kable_styling(latex_options=c("hold_position", font_size=11)) %>%\n  add_footnote(c("information from samtools flagstat"))\n\nreadCountTargets <- file.path("Analysis", "Minimap", \n                              paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE), ".bam", sep=""))\n\nExternalAnnotation = file.path("ReferenceData", basename(config$genome_annotation))\n\ngeneCounts <- featureCounts(files=readCountTargets,\n                            annot.ext=ExternalAnnotation,\n                            isGTFAnnotationFile=TRUE,\n                            GTF.featureType="exon",\n                            GTF.attrType="gene_id",\n                            isLongRead=TRUE,\n                            largestOverlap=TRUE,\n                            useMetaFeatures=TRUE,\n                            nthreads = 8)$counts\n\ncolnames(geneCounts) <- rownames(studyDesign)\n\nknitr::kable(geneCounts[order(rowSums(geneCounts), decreasing=TRUE)[1:10],], caption="Table showing the 10 annotated gene features with the highest number of mapped reads", booktabs=TRUE, table.envir=\'table*\', linesep="") %>%\n  kable_styling(latex_options=c("hold_position", font_size=11)) %>%\n  add_footnote(c("This is raw count data and no normalisation or transformation has been performed"))\n\ngeneCounts_nonZeros <- geneCounts[which(rowSums(geneCounts) > 0),] \n\nens.geneCounts_nonZeros <- rownames(geneCounts_nonZeros)\nSymbols <- mapIds(org.Hs.eg.db, keys = ens.geneCounts_nonZeros, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")\n\ngeneCounts_nonZeros <- as.data.frame(geneCounts_nonZeros)\ngeneCounts_nonZeros$Symbols <- Symbols\n\nxlsExpressedGenes <- file.path(resultDir, "ExpressedGenes.xlsx")\ngeneCounts_nonZeros$gene_id <- rownames(geneCounts_nonZeros) \nwrite_xlsx(x = geneCounts_nonZeros, path = xlsExpressedGenes)\n\ndeSeqRaw <- DESeqDataSetFromMatrix(countData=geneCounts, colData=studyDesign, design= ~group)\n\ndeSeqRaw <- deSeqRaw[-which(rowSums(counts(deSeqRaw)) < config$readCountMinThreshold), ]\n\ncolnames(geneCounts) <- studyDesign$group\n\ngeneCounts_bar <- geneCounts\ndim(geneCounts_bar)\n\nens.geneCounts_bar <- rownames(geneCounts_bar)\nSymbols_bar <- mapIds(org.Hs.eg.db, keys = ens.geneCounts_bar, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")\ngeneCounts_bar <- as.data.frame(geneCounts_bar)\ngeneCounts_bar$Symbols <- Symbols_bar\n\nis.na(geneCounts_bar$Symbols)\ngeneCounts_bar <- drop_na(geneCounts_bar)\ndim(geneCounts_bar)\n\ngeneCounts_bar <- geneCounts_bar[!duplicated(geneCounts_bar$Symbols),]\n\nrownames(geneCounts_bar) <- geneCounts_bar$Symbols\ngeneCounts_bar$Symbols <- NULL\ngeneCounts_bar <- as.data.frame(t(geneCounts_bar))\n\ngeneCounts_bar$Treatment <- rownames(geneCounts_bar)\ngeneCounts_bar$Treatment <- factor(geneCounts_bar$Treatment, levels = c("barcode05", "barcode06", "barcode07", "barcode08",  "barcode09", "barcode11"))\n\nggplot(geneCounts_bar)+\n  geom_col(aes(x = Treatment, y = AKR1C3, fill = Treatment))+\n  ggtitle("Raw counts: AKR1C3")+\n  scale_fill_brewer(palette = "Paired")+\n  theme(axis.text.x = element_text(angle = 90, hjust = 1))'







def group_set_barcodes(sq_list, s_id, barcode_fw = 'static/test/barcode_fw.fasta', barcode_rv = 'static/test/barcode_rw.fasta'):
    
    barcodes_set_fw = list(SeqIO.parse(barcode_fw, 'fasta'))
    barcodes_set_rv = list(SeqIO.parse(barcode_rv, 'fasta'))
    
    # !TODO:PARALLEL THIS PART
    for barcode_fw in barcodes_set_fw:
        for barcode_rv in barcodes_set_rv:
            idx1 = barcode_fw.description
            idx2 = barcode_rv.description
            groupby_barcodes(sq_list, s_id, idx= idx1 +'_' + idx2, barcode_fw= barcode_fw.seq, barcode_rv = barcode_rv.seq)
            run_minimap2(s_id=s_id, idx= idx1 +'_' + idx2)
    shutil.make_archive('static/%s'%s_id, 'zip', 'static/test_result/%s' %s_id)
    return

    # May remove upload data in 30 days


# 
def send_result(submission_name, link_address, recipient_email):
    send_mail(
        'Your submission < > completed',
        'Your computational result is stored within 30 days, You can download your result from this link',
        'thach.iuf@gmail.com',
        recipient_list=[recipient_email],
        fail_silently=False,
    )
    return

# Start WT and genome alignment
# Require pyensembl, minimap2
# pyensembl install --release 102 --species homo_sapiens




def get_genome_data(gene_name = 'ACE2'):
    from pyensembl import EnsemblRelease
    release = EnsemblRelease()[0]


if __name__ == "__main__":
    import doctest
    doctest.testmod()




'setwd("/home/ag-rossi/projects/NGS_webserver/NGS_webserver/RNAseq_analyzer/%s")\n\
    library(digest)\n\
    library(ShortRead)\n\
    library(ggplot2)\n\
    library(plyr)\n\
    library(dplyr)\n\
    library(tidyr)\n\
    library(Rsubread)\n\
    library(DESeq2)\n\
    library(pcaMethods)\n\
    library(kableExtra)\n\
    library(caTools)\n\
    library(writexl)\n\
    library(yaml)\n\
    library(session)\n\
    library(AnnotationDbi)\n\
    library(org.Hs.eg.db)\n\
    source("Static/R/common.R")\n\
    resultDir <- file.path("Analysis", "Results")\n\
    dir.create(resultDir, showWarnings = FALSE, recursive = TRUE)\n\
    config <- yaml::yaml.load_file("config.yaml")\n\
    persistenceData <- file.path(resultDir, "NanoporeDESeq2.Rdata")\n\
    studyDesign <- data.frame() \nfor (i in 1:length(config$Samples)) {\n  studyDesign <- rbind(studyDesign,\n                       data.frame(samples=names(config$Samples[[i]][[1]]), \n                                  filename=unlist(config$Samples[[i]][[1]]), \n                                  group=names(config$Samples[[i]])))\n}\n\
    studyDesign$replicate <- sapply(1:nrow(studyDesign), function(x)sum(studyDesign$group[1:x]==studyDesign$group[x]))\n\
    studyDesign$md5 <- lapply(as.character(studyDesign$filename), md5sum)\n\
    studyDesign <- studyDesign[,-which(colnames(studyDesign)=="samples")]\n\
    processQCFastq <- function(rowname) {\n  row <- which(row.names(studyDesign)==rowname)\n  file <- as.character(studyDesign[row, "filename"])\n  fastq <- readFastq(file, widthIds=FALSE)\n  c(\n    reads = formatC(length(fastq), big.mark=","),\n    mbs = formatC(round(sum(width(fastq)) / 1000 / 1000, digits=1), big.mark=","),\n    min = min(width(fastq)),\n    max = max(width(fastq)),\n    mean = round(mean(width(fastq)), digits=1),\n    median = round(median(width(fastq)), digits=0),\n    qval = round(mean(alphabetScore(fastq) / width(fastq)), digits=1),\n    gc = round(mean(letterFrequency(sread(fastq), "GC")  / width(fastq)) * 100, digits=1),\n    n50 = ncalc(width(fastq), n=0.5), # The ncalc function is defined in the common.R-File\n    l50 = lcalc(width(fastq), n=0.5),\n    n90 = ncalc(width(fastq), n=0.9),\n    l90 = lcalc(width(fastq), n=0.9)\n  )\n}\n\
    data <- lapply(row.names(studyDesign), processQCFastq)\n\
    qcData <- data.frame(data)\n\
    colnames(qcData) <- row.names(studyDesign)\n\
    View(qcData)\n\
    knitr::kable(qcData, caption="Summary statistics for the cDNA libraries imported", booktabs=TRUE, table.envir=\'table*\', linesep="")  %>%\n  kable_styling(latex_options=c("hold_position", font_size=9))\n\
    extractLengths <- function(rowname) {\n  row <- which(row.names(studyDesign)==rowname)\n  file <- as.character(studyDesign[row, "filename"])\n  fastq <- readFastq(file)\n  t(as.matrix(width(fastq)))\n}\n\
    lengthData <- mapply(extractLengths, row.names(studyDesign))\n\
    lengthDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(lengthData)))\ncolnames(lengthDataMatrix) <-  row.names(studyDesign)\n\n\
    lengthMatrixMelt <- reshape2::melt(lengthDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))\n\
    lengthMatrixMelt <- cbind(lengthMatrixMelt, group=studyDesign[match(lengthMatrixMelt$variable,  rownames(studyDesign)), "group"])\n\n\
    plot <- ggplot(lengthMatrixMelt, aes(x=variable, y=value, fill=group)) + geom_violin() + scale_y_continuous(limits=c(0, as.numeric(quantile(lengthMatrixMelt$value, probs=c(0.975))))) + xlab("study sample") +  ylab("Distribution of Read Lengths (bp)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution of read lengths across samples")\n\n\
    suppressWarnings(print(plot))\n\n\
    extractQualities <- function(rowname) {\n  row <- which(row.names(studyDesign)==rowname)\n  file <- as.character(studyDesign[row, "filename"])\n  fastq <- readFastq(file)\n  t(as.matrix(alphabetScore(fastq) / width(fastq)))\n}\n\
    qualityData <- mapply(extractQualities, row.names(studyDesign))\n\
    qualityDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(qualityData)))\n\
    colnames(qualityDataMatrix) <-  row.names(studyDesign)\n\n\
    qualityMatrixMelt <- reshape2::melt(qualityDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))\n\
    qualityMatrixMelt <- cbind(qualityMatrixMelt, group=studyDesign[match(qualityMatrixMelt$variable,  rownames(studyDesign)), "group"])\n\n\
    plotQ <- ggplot(qualityMatrixMelt, aes(x=variable, y=value, fill=group)) + geom_violin() + scale_y_continuous(limits=c(min(qualityMatrixMelt$value), max(qualityMatrixMelt$value))) + xlab("study sample") +  ylab("Distribution of Read Qualities (QV)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution of read qualities across samples")\n\n\
    suppressWarnings(print(plotQ))\n\n\
    flagstatTargets <- file.path("Analysis", "flagstat", \n                             paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE),".txt",sep=""))\n\n\
    loadFlagstat <- function(file) {\n  x <- read.table(file, header=FALSE, sep=" ", fill=NA)[c(1:5),c(1,3)]\n  x[,1]\n}\n\n\
    flagstatRes <- data.frame(matrix(unlist(lapply(flagstatTargets, loadFlagstat)), ncol=length(flagstatTargets)), stringsAsFactors = FALSE)\n\
    colnames(flagstatRes) <- rownames(studyDesign)\n\
    rownames(flagstatRes) <- c("read mappings", "Secondary", "Supplementary", "Duplicates", "Mapped")\n\n\
    flagstatRes[nrow(flagstatRes)+1,] <- as.numeric(gsub(",","",t(qcData)[, "reads"]))\n\
    rownames(flagstatRes)[6] <- "nreads"\n\n\
    getVal <- function(word) {\n  sum(as.numeric(unlist(strsplit(word, "/"))))\n}\n\n\
    zreads <- unlist(lapply(flagstatRes["read mappings", ], getVal)) -\n  unlist(lapply(flagstatRes["Secondary", ], getVal)) -\n  unlist(lapply(flagstatRes["Supplementary", ], getVal)) - \n  unlist(lapply(flagstatRes["Duplicates", ], getVal)) \n\n\
    flagstatRes[nrow(flagstatRes)+1,] <- round(as.numeric(flagstatRes["Mapped", ]) / as.numeric(flagstatRes["nreads", ]) * 100, digits = 2)\n\n\
    flagstatRes <- flagstatRes[c(6,1,2,3,4,7),]\n\
    rownames(flagstatRes)[6] <- "%mapping"\n\nknitr::kable(flagstatRes, caption="Summary statistics from the minimap2 long read spliced mapping.", booktabs=TRUE, table.envir=\'table*\', linesep="")  %>%\n  kable_styling(latex_options=c("hold_position", font_size=11)) %>%\n  add_footnote(c("information from samtools flagstat"))\n\nreadCountTargets <- file.path("Analysis", "Minimap", \n                              paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE), ".bam", sep=""))\n\nExternalAnnotation = file.path("ReferenceData", basename(config$genome_annotation))\n\n\
    geneCounts <- featureCounts(files=readCountTargets,\n                            annot.ext=ExternalAnnotation,\n                            isGTFAnnotationFile=TRUE,\n                            GTF.featureType="exon",\n                            GTF.attrType="gene_id",\n                            isLongRead=TRUE,\n                            largestOverlap=TRUE,\n                            useMetaFeatures=TRUE,\n                            nthreads = 8)$counts\n\n\
    colnames(geneCounts) <- rownames(studyDesign)\n\n\
    knitr::kable(geneCounts[order(rowSums(geneCounts), decreasing=TRUE)[1:10],], caption="Table showing the 10 annotated gene features with the highest number of mapped reads", booktabs=TRUE, table.envir=\'table*\', linesep="") %>%\n  kable_styling(latex_options=c("hold_position", font_size=11)) %>%\n  add_footnote(c("This is raw count data and no normalisation or transformation has been performed"))\n\ngeneCounts_nonZeros <- geneCounts[which(rowSums(geneCounts) > 0),] \n\nens.geneCounts_nonZeros <- rownames(geneCounts_nonZeros)\nSymbols <- mapIds(org.Hs.eg.db, keys = ens.geneCounts_nonZeros, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")\n\ngeneCounts_nonZeros <- as.data.frame(geneCounts_nonZeros)\ngeneCounts_nonZeros$Symbols <- Symbols\n\nxlsExpressedGenes <- file.path(resultDir, "ExpressedGenes.xlsx")\ngeneCounts_nonZeros$gene_id <- rownames(geneCounts_nonZeros) \nwrite_xlsx(x = geneCounts_nonZeros, path = xlsExpressedGenes)\n\ndeSeqRaw <- DESeqDataSetFromMatrix(countData=geneCounts, colData=studyDesign, design= ~group)\n\ndeSeqRaw <- deSeqRaw[-which(rowSums(counts(deSeqRaw)) < config$readCountMinThreshold), ]\n\ncolnames(geneCounts) <- studyDesign$group\n\ngeneCounts_bar <- geneCounts\ndim(geneCounts_bar)\n\nens.geneCounts_bar <- rownames(geneCounts_bar)\nSymbols_bar <- mapIds(org.Hs.eg.db, keys = ens.geneCounts_bar, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")\ngeneCounts_bar <- as.data.frame(geneCounts_bar)\ngeneCounts_bar$Symbols <- Symbols_bar\n\nis.na(geneCounts_bar$Symbols)\ngeneCounts_bar <- drop_na(geneCounts_bar)\ndim(geneCounts_bar)\n\ngeneCounts_bar <- geneCounts_bar[!duplicated(geneCounts_bar$Symbols),]\n\nrownames(geneCounts_bar) <- geneCounts_bar$Symbols\ngeneCounts_bar$Symbols <- NULL\ngeneCounts_bar <- as.data.frame(t(geneCounts_bar))\n\ngeneCounts_bar$Treatment <- rownames(geneCounts_bar)\ngeneCounts_bar$Treatment <- factor(geneCounts_bar$Treatment, levels = c("barcode05", "barcode06", "barcode07", "barcode08",  "barcode09", "barcode11"))\n\n\
    ggplot(geneCounts_bar)+\n  geom_col(aes(x = Treatment, y = AKR1C3, fill = Treatment))+\n  ggtitle("Raw counts: AKR1C3")+\n  scale_fill_brewer(palette = "Paired")+\n  theme(axis.text.x = element_text(angle = 90, hjust = 1))'%(path)
