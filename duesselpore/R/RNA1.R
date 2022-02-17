# Loading required packages
suppressMessages(library(digest))
suppressMessages(library(ShortRead))  
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))      
suppressMessages(library(tidyr))
suppressMessages(library(Rsubread))
suppressMessages(library(DESeq2))
suppressMessages(library(pcaMethods))
suppressMessages(library(kableExtra))  
suppressMessages(library(caTools))     
suppressMessages(library(writexl))     
suppressMessages(library(yaml))
suppressMessages(library(session))
suppressMessages(library(AnnotationDbi)) 
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(RColorBrewer))
suppressMessages(library(PoiClaClu))
suppressMessages(library(pheatmap))
suppressMessages(library(karyoploteR))
suppressMessages(library(ensembldb))
suppressMessages(library(viridis))
suppressMessages(library(gridExtra))
suppressMessages(library(DOSE))
suppressMessages(library(pathview))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))

suppressMessages(library(EnsDb.Hsapiens.v86))
refdb<-EnsDb.Hsapiens.v86

# Custamized functions----
source("Static/R/common.R")

resultDir <- file.path("Analysis", "Results")
dir.create(resultDir, showWarnings = FALSE, recursive = TRUE) 

# Loading the config file----
config <- yaml::yaml.load_file("config.yaml")

# Defining the file where the session and other data will be saved----
persistenceData <- file.path(resultDir, "NanoporeDESeq2.Rdata")

# Create a study design by using the config.yaml-file
studyDesign <- data.frame()
for (i in 1:length(config$Samples)) {
  studyDesign <- rbind(studyDesign,
                       data.frame(samples=names(config$Samples[[i]][[1]]), 
                                  filename=unlist(config$Samples[[i]][[1]]), 
                                  group=names(config$Samples[[i]])))
}

studyDesign$replicate <- sapply(1:nrow(studyDesign), function(x)sum(studyDesign$group[1:x]==studyDesign$group[x]))
studyDesign$md5 <- lapply(as.character(studyDesign$filename), md5sum)
#studyDesign <- studyDesign[,-which(colnames(studyDesign)=="samples")]

# Raw sequence review----
processQCFastq <- function(rowname) {
  row <- which(row.names(studyDesign)==rowname)
  file <- as.character(studyDesign[row, "filename"])
  fastq <- readFastq(file, widthIds=FALSE)
  c(
    reads = formatC(length(fastq), big.mark=","),
    mbs = formatC(round(sum(width(fastq)) / 1000 / 1000, digits=1), big.mark=","),
    min = min(width(fastq)),
    max = max(width(fastq)),
    mean = round(mean(width(fastq)), digits=1),
    median = round(median(width(fastq)), digits=0),
    qval = round(mean(alphabetScore(fastq) / width(fastq)), digits=1),
    gc = round(mean(letterFrequency(sread(fastq), "GC")  / width(fastq)) * 100, digits=1)
    #n50 = ncalc(width(fastq), n=0.5),
    #l50 = lcalc(width(fastq), n=0.5),
    #n90 = ncalc(width(fastq), n=0.9),
    #l90 = lcalc(width(fastq), n=0.9)
  )
}

data <- lapply(row.names(studyDesign), processQCFastq)
qcData <- data.frame(data)
colnames(qcData) <- row.names(studyDesign)
pdf("Analysis/Results/sample_summary.pdf", width=12, height = 4)       # Export PDF
grid.table(qcData)
dev.off()

# Distribution of read Lengths (bp)----
extractLengths <- function(rowname) {
  row <- which(row.names(studyDesign)==rowname)
  file <- as.character(studyDesign[row, "filename"])
  fastq <- readFastq(file)
  t(as.matrix(width(fastq)))
}
lengthData <- mapply(extractLengths, row.names(studyDesign))
lengthDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(lengthData)))
colnames(lengthDataMatrix) <-  row.names(studyDesign)

lengthMatrixMelt <- reshape2::melt(lengthDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))
lengthMatrixMelt <- cbind(lengthMatrixMelt, group=studyDesign[match(lengthMatrixMelt$variable,  rownames(studyDesign)), "group"])

plot <- ggplot(lengthMatrixMelt, aes(x=variable, y=value, fill=group)) + geom_violin() + scale_y_continuous(limits=c(0, as.numeric(quantile(lengthMatrixMelt$value, probs=c(0.975))))) + xlab("study sample") +  ylab("Distribution of Read Lengths (bp)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution of read lengths across samples")
ggsave("Analysis/Results/read_length.pdf", width=12, height = 6)

# read quality----
extractQualities <- function(rowname) {
  row <- which(row.names(studyDesign)==rowname)
  file <- as.character(studyDesign[row, "filename"])
  fastq <- readFastq(file)
  t(as.matrix(alphabetScore(fastq) / width(fastq)))
}
qualityData <- mapply(extractQualities, row.names(studyDesign))
qualityDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(qualityData)))
colnames(qualityDataMatrix) <-  row.names(studyDesign)

qualityMatrixMelt <- reshape2::melt(qualityDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))
qualityMatrixMelt <- cbind(qualityMatrixMelt, group=studyDesign[match(qualityMatrixMelt$variable,  rownames(studyDesign)), "group"])

plotQ <- ggplot(qualityMatrixMelt, aes(x=variable, y=value, fill=group)) + geom_violin() + scale_y_continuous(limits=c(min(qualityMatrixMelt$value), max(qualityMatrixMelt$value))) + xlab("study sample") +  ylab("Distribution of Read Qualities (QV)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution of read qualities across samples")
ggsave("Analysis/Results/read_quality.pdf", width=12, height = 6)


# Review of cDNA read mapping----
# Read mappings correspond to the number of unique reads. In our worklflow we disabled secondary mappings.
# A Supplementary alignment could correspond to an chimeric read. For instance a result of a read with a barcode located in the middle.
# The informations are obtained from the .txt-Files with the flagstats, created from the .bam-Files.
flagstatTargets <- file.path("Analysis", "flagstat", 
                             paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE),".txt",sep=""))

loadFlagstat <- function(file) {
  x <- read.table(file, header=FALSE, sep=" ", fill=NA)[c(1:5),c(1,3)]
  x[,1]
}

flagstatRes <- data.frame(matrix(unlist(lapply(flagstatTargets, loadFlagstat)), ncol=length(flagstatTargets)), stringsAsFactors = FALSE)
colnames(flagstatRes) <- rownames(studyDesign)
rownames(flagstatRes) <- c("read mappings", "Secondary", "Supplementary", "Duplicates", "Mapped")


flagstatRes[nrow(flagstatRes)+1,] <- as.numeric(gsub(",","",t(qcData)[, "reads"]))
rownames(flagstatRes)[6] <- "nreads"

getVal <- function(word) {
  sum(as.numeric(unlist(strsplit(word, "/"))))
}

zreads <- unlist(lapply(flagstatRes["read mappings", ], getVal)) -
  unlist(lapply(flagstatRes["Secondary", ], getVal)) -
  unlist(lapply(flagstatRes["Supplementary", ], getVal)) - 
  unlist(lapply(flagstatRes["Duplicates", ], getVal)) 

flagstatRes[nrow(flagstatRes)+1,] <- round(as.numeric(flagstatRes["Mapped", ]) / as.numeric(flagstatRes["nreads", ]) * 100, digits = 2)

flagstatRes <- flagstatRes[c(6,1,2,3,4,7),]
rownames(flagstatRes)[6] <- "%mapping"

pdf("Analysis/Results/cDNA_mapping_summary.pdf", width=12, height = 3)       # Export PDF
grid.table(flagstatRes)
dev.off()
