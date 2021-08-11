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
Symbols <- mapIds(refdb, keys = ens.geneCounts_nonZeros, column = "SYMBOL", keytype = "GENEID", multiVals = "first")

# To add the Symbols in an additional column, the matrix needs to be converted into a DataFrame
geneCounts_nonZeros <- as.data.frame(geneCounts_nonZeros)
geneCounts_nonZeros$Symbols <- Symbols

# Save the raw count data in an Excel-file
xlsExpressedGenes <- file.path(resultDir, "ExpressedGenes.xlsx")
geneCounts_nonZeros$gene_id <- rownames(geneCounts_nonZeros) # As Rownames will not be exported, we create an additional column 
write_xlsx(x = geneCounts_nonZeros, path = xlsExpressedGenes)


