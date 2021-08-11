library("readxl")
count_data <-read_excel('Analysis/Results/ExpressedGenes.xlsx')
geneCounts <-count_data %>% select(studyDesign$samples)
colnames(geneCounts) <-studyDesign$samples

geneCounts <-as.matrix(geneCounts) 
row.names(geneCounts)<-c(count_data$gene_id)
geneCounts_nonZeros <- geneCounts[which(rowSums(geneCounts) > 0),] 
# Using the AnnotationID package to convert the ENSEMBLE ID into Gene symbols. Therefore we use the org.Hs.eg.dg Database
ens.geneCounts_nonZeros <- rownames(geneCounts_nonZeros)
Symbols <- mapIds(refdb, keys = ens.geneCounts_nonZeros, column = "SYMBOL", keytype = "GENEID", multiVals = "first")

# To add the Symbols in an additional column, the matrix needs to be converted into a DataFrame
geneCounts_nonZeros <- as.data.frame(geneCounts_nonZeros)
geneCounts_nonZeros$Symbols <- Symbols


xlsExpressedGenes <- file.path(resultDir, "ExpressedGenes.xlsx")
geneCounts_nonZeros$gene_id <- rownames(geneCounts_nonZeros) # As Rownames will not be exported, we create an additional column 
write_xlsx(x = geneCounts_nonZeros, path = xlsExpressedGenes)



