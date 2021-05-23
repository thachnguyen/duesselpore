library(digest)
library(ShortRead)  
library(ggplot2)
library(plyr)
library(dplyr)      
library(tidyr)
library(Rsubread)
library(DESeq2)
library(pcaMethods)
library(kableExtra)  
library(caTools)     # base64decode for embedded Snakefile
library(writexl)     # writing results to Excel files
library(yaml)
library(session)
library(AnnotationDbi) # These packages are needed to convert the ENSEMBL ID into Gene symbols
library(org.Hs.eg.db)

# Loading of a file with customized functions----
source("Static/R/common.R")

# Creating of a folder in which the results will be stored----
resultDir <- file.path("Analysis", "Results")
dir.create(resultDir, showWarnings = FALSE, recursive = TRUE) #recursive defines whether other folder of this path are supposd to be created

# Defining the configuration----
# The configuration can be predefined in a .yaml-File
# To use this, it has to be loaded before the analysis
config <- yaml::yaml.load_file("config.yaml")

# Defining the file where the session and other data will be saved----
persistenceData <- file.path(resultDir, "NanoporeDESeq2.Rdata")

# Create a study design by using the config.yaml-File
studyDesign <- data.frame() # hierfür benötigen wir zunächst ein leeres DataFrame
#To fill this DataFrame with some values, we create a for loop, which iterates through the config-file to obtain the data
for (i in 1:length(config$Samples)) {
  studyDesign <- rbind(studyDesign,
                       data.frame(samples=names(config$Samples[[i]][[1]]), 
                                  filename=unlist(config$Samples[[i]][[1]]), 
                                  group=names(config$Samples[[i]])))
}
# With this function we add a new columns containing informations about the number of replicates
studyDesign$replicate <- sapply(1:nrow(studyDesign), function(x)sum(studyDesign$group[1:x]==studyDesign$group[x]))

# Adding a md5 fingerprint.
# Theoretically there are unlimited combinations and number of files with any given MD5-Hash
# Basically a control that the accidentally aligned a fastq-File twise. This would result in a identical MD5-Hash
studyDesign$md5 <- lapply(as.character(studyDesign$filename), md5sum)

# Relevel the data according to our control group.
#studyDesign$group <- relevel(studyDesign$group, ref=config$referenceGroup)

# Removing the column "Samples" as it became redundant. Tidy data!!
studyDesign <- studyDesign[,-which(colnames(studyDesign)=="samples")]

# Alternatively a txt-File can be prepared in advanced to load in the study desing, but we want to automate as many steps as possible

# Raw sequence review----
# In this step we try to get a deeper view into the raw sequencing data.
# For this we make use of the ShortRead-Package
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
    gc = round(mean(letterFrequency(sread(fastq), "GC")  / width(fastq)) * 100, digits=1),
    n50 = ncalc(width(fastq), n=0.5), # The ncalc function is defined in the common.R-File
    l50 = lcalc(width(fastq), n=0.5),
    n90 = ncalc(width(fastq), n=0.9),
    l90 = lcalc(width(fastq), n=0.9)
  )
}

data <- lapply(row.names(studyDesign), processQCFastq)
qcData <- data.frame(data)
colnames(qcData) <- row.names(studyDesign)
View(qcData)

# Using the package knitr to create a noce and dynamic report
knitr::kable(qcData, caption="Summary statistics for the cDNA libraries imported", booktabs=TRUE, table.envir='table*', linesep="")  %>%
  kable_styling(latex_options=c("hold_position", font_size=9))
# We can use the number of reads, overall read quality and GC content to get an impression about the differences of the samples. 
# Huge differences may indicate technical differences and might be problematic for DE-analysis

# Using a violin plot to show the distribution of read Lengths (bp)----
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

suppressWarnings(print(plot))

# Violin plot show the quality of reads----
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

suppressWarnings(print(plotQ))

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

# Read mappings are the product of all alignments. In this case the number of reads plus the Supplementary alignments.
# The number of reads will also be included. These are the real read numbers.
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

#knitr::kable(flagstatRes, caption="Summary statistics from the Minimap2 long read spliced mapping.")

knitr::kable(flagstatRes, caption="Summary statistics from the minimap2 long read spliced mapping.", booktabs=TRUE, table.envir='table*', linesep="")  %>%
  kable_styling(latex_options=c("hold_position", font_size=11)) %>%
  add_footnote(c("information from samtools flagstat"))

# Analysis of reads mapped to genes----
# Here we use the featureCount method from the package Rsubread
# For this the .bam-Files are required. In a first step the path to this files will be defined.
readCountTargets <- file.path("Analysis", "Minimap", 
                              paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE), ".bam", sep=""))

# Rsubread already has some inbuilt Annotations. But in this case we define a external Annotation, which was used previously.
ExternalAnnotation = file.path("/home/ag-rossi/ReferenceData", basename(config$genome_annotation))

# With the Method featureCounts, the reads will be mapped to genes
# This step might take a while, Get a coffee!!
geneCounts <- featureCounts(files=readCountTargets,
                            annot.ext=ExternalAnnotation,
                            isGTFAnnotationFile=TRUE,
                            GTF.featureType="exon",
                            GTF.attrType="gene_id",
                            isLongRead=TRUE,
                            largestOverlap=TRUE,
                            useMetaFeatures=TRUE,
                            nthreads = 16)$counts
# largestOverlap sorts the multi-mapping challenge with long reads
# reportReads="CORE" can be appended to get info on read status

# Rename column headers for clarity
colnames(geneCounts) <- rownames(studyDesign)

# Creating a table with genes showing the highest number of counts----
# Keep in mind that neither a normalization nor a transformation was performed 
knitr::kable(geneCounts[order(rowSums(geneCounts), decreasing=TRUE)[1:10],], caption="Table showing the 10 annotated gene features with the highest number of mapped reads", booktabs=TRUE, table.envir='table*', linesep="") %>%
  kable_styling(latex_options=c("hold_position", font_size=11)) %>%
  add_footnote(c("This is raw count data and no normalisation or transformation has been performed"))

# Convert the ENSEMBL ID into Gene symbols----
# Removing those genes without a single read
geneCounts_nonZeros <- geneCounts[which(rowSums(geneCounts) > 0),] 

# Using the AnnotationID package to convert the ENSEMBLE ID into Gene symbols. Therefore we use the org.Hs.eg.dg Database
ens.geneCounts_nonZeros <- rownames(geneCounts_nonZeros)
Symbols <- mapIds(org.Hs.eg.db, keys = ens.geneCounts_nonZeros, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# To add the Symbols in an additional column, the matrix needs to be converted into a DataFrame
geneCounts_nonZeros <- as.data.frame(geneCounts_nonZeros)
geneCounts_nonZeros$Symbols <- Symbols

# Save the raw count data in an Excel-file
xlsExpressedGenes <- file.path(resultDir, "ExpressedGenes.xlsx")
geneCounts_nonZeros$gene_id <- rownames(geneCounts_nonZeros) # As Rownames will not be exported, we create an additional column 
write_xlsx(x = geneCounts_nonZeros, path = xlsExpressedGenes)

# Analyzing the differential expression using DESeq2
# DESeq2 utilises the raw read count data to model between condition variability and to identify the differentially expressed genes. 
# DESeq2 does not use normalized data, but corrects internally for the relative library size to assess measurement precision. 
# Thresholds have been defined in the config.yaml-File

deSeqRaw <- DESeqDataSetFromMatrix(countData=geneCounts, colData=studyDesign, design= ~group)

keep <- rowSums(counts(deSeqRaw)) > config$readCountMinThreshold
deSeqRaw <- deSeqRaw[keep,]
dds <- DESeq(deSeqRaw)
# filter out the features that do not contain at least min threshold of expressed genes


vsd <- vst(object = deSeqRaw, blind = TRUE)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), n = 50)
mat <- assay(vsd)[topVarGenes,]
mat <- mat - rowMeans(mat)

mat_breaks <- seq(min(mat), max(mat), length.out = 10)
anno <- as.data.frame(colData(vsd)[, c("group","replicate")])
anno$replicate <- factor(anno$replicate) 
variance_heatmap <- pheatmap(mat, 
         annotation_col = anno, 
         labels_row = mapIds(org.Hs.eg.db, keys = substr(rownames(mat),1,15), 
                             column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
         labels_col = c(studyDesign$group, studyDesign$replicate), drop_levels = TRUE, filename='heatmap.pdf')

organism <- org.Hs.eg.db
options(repr.plot.width=12, repr.plot.height=15)



gse <- gseGO(geneList = geneList,
             ont = "ALL",
             keyType = "ENSEMBL",
             #nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "none",
             eps =0.0,
             filename = 'pathway.pdf'
            )

dotplot(gse, showCategory = 35, split =".sign", orderBy = "x")+
  facet_grid(.~.sign)


