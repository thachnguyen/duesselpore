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
studyDesign <- studyDesign[,-which(colnames(studyDesign)=="samples")]

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
    gc = round(mean(letterFrequency(sread(fastq), "GC")  / width(fastq)) * 100, digits=1),
    n50 = ncalc(width(fastq), n=0.5),
    l50 = lcalc(width(fastq), n=0.5),
    n90 = ncalc(width(fastq), n=0.9),
    l90 = lcalc(width(fastq), n=0.9)
  )
}

data <- lapply(row.names(studyDesign), processQCFastq)
qcData <- data.frame(data)
colnames(qcData) <- row.names(studyDesign)
View(qcData)

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
ggsave("Analysis/Results/read_length.pdf")

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
ggsave("Analysis/Results/read_quality.pdf")


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


# Analysis of reads mapped to genes----
readCountTargets <- file.path("Analysis", "Minimap", 
                              paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE), ".bam", sep=""))

# Rsubread already has some inbuilt Annotations. But in this case we define a external Annotation, which was used previously.
ExternalAnnotation = file.path("/home/ag-rossi/ReferenceData", basename(config$genome_annotation))

# With the Method featureCounts, the reads will be mapped to genes
geneCounts <- featureCounts(files=readCountTargets,
                            annot.ext=ExternalAnnotation,
                            isGTFAnnotationFile=TRUE,
                            GTF.featureType="exon",
                            GTF.attrType="gene_id",
                            isLongRead=TRUE,
                            largestOverlap=TRUE,
                            useMetaFeatures=TRUE,
                            nthreads = 16)$counts

# Rename column headers for clarity
colnames(geneCounts) <- rownames(studyDesign)

# Creating a table with genes showing the highest number of counts----
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
group_size <- 2
for (group1 in unique(studyDesign$group)) {
  if (group_size > sum(studyDesign$group== group1))
      {group_size <- 1}
  }
if (group_size > 1) { 
  deSeqRaw <- DESeqDataSetFromMatrix(countData=geneCounts, colData=studyDesign, design= ~group)
} else {
  deSeqRaw <- DESeqDataSetFromMatrix(countData=geneCounts, colData=studyDesign, design= ~ 1)  
  }
keep <- rowSums(counts(deSeqRaw)) > config$readCountMinThreshold
deSeqRaw <- deSeqRaw[keep,]
dds <- DESeq(deSeqRaw)
# filter out the features that do not contain at least min threshold of expressed genes

if (config$organism=='human'){
vsd <- vst(object = deSeqRaw, blind = TRUE)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), n = 50)
mat <- assay(vsd)[topVarGenes,]
mat <- mat - rowMeans(mat)

mat_breaks <- seq(min(mat), max(mat), length.out = 10)
anno <- as.data.frame(colData(vsd)[, c("group","replicate")])
anno$replicate <- factor(anno$replicate) 
variance_heatmap <- pheatmap(mat, 
         annotation_col = anno, 
         cluster_cols=config$cluster_col,
         labels_row = mapIds(org.Hs.eg.db, keys = substr(rownames(mat),1,15), 
                             column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
         labels_col = c(studyDesign$group, studyDesign$replicate), drop_levels = TRUE, filename='Analysis/Results/heatmap.pdf')

organism <- org.Hs.eg.db
res_group01_group02.filtered <- res_group01_group02 %>%
  as.data.frame() %>%
  dplyr::filter(abs(log2FoldChange) > 1.5)

geneList <- res_group01_group02.filtered$log2FoldChange
names(geneList) <- rownames(res_group01_group02.filtered)
geneList <- sort(geneList, decreasing = TRUE)

options(repr.plot.width=12, repr.plot.height=15)

gse <- gseGO(geneList = geneList,
             ont = "ALL",
             keyType = "ENSEMBL",
             #nPerm = 10000,
             #minGSSize = 5,
             #maxGSSize = 500,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "none",
             #eps =0.0,
)

dotplot(gse, showCategory = 35, split =".sign", orderBy = "x")+
  facet_grid(.~.sign)
ggsave("Analysis/Results/gseGO.pdf", width=20, height = 15)

gse.df <- as.data.frame(gse)
pathway <- file.path("Analysis/Results/pathway.xlsx")
write_xlsx(x= gse.df, path = pathway)

data(geneList)
de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)
library(enrichplot)
barplot(edo, showCategory=20)
ggsave("Analysis/Results/bar_plot_enrichDGN.pdf")

options(repr.plot.width=20, repr.plot.height=7)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)

## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
ggsave("Analysis/Results/cnet_plot1.pdf")

options(repr.plot.width=20, repr.plot.height=12)
p1 <- cnetplot(edox, node_label="category") 
p2 <- cnetplot(edox, node_label="gene") 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
ggsave("Analysis/Results/cnet_plot2.pdf")

}


