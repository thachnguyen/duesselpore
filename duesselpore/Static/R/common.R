
md5sum <- function(filename) digest::digest(filename, algo="md5", file=TRUE)

ncalc <- function(len.vector, n) {
  # N50 - length such that scaffolds of this length or longer include half the bases of the assembly
  len.sorted <- rev(sort(len.vector))
  len.sorted[cumsum(len.sorted) >= sum(len.sorted)*n][1]
}


lcalc <- function(len.vector, n) {
  len.sorted <- rev(sort(len.vector))
  which(cumsum(len.sorted) >= sum(len.sorted)*n)[1]
}


volcanoPlot <- function(deSeqObj, lfcThreshold, adjPValueThreshold) {
  volcanoSubstr <- results(deSeqObj)
  logUp <- which(volcanoSubstr$log2FoldChange >= lfcThreshold)
  logDown <- which(volcanoSubstr$log2FoldChange <= -lfcThreshold)
  withStat <- which(volcanoSubstr$padj <= adjPValueThreshold)
  colours <- c(noDifference="gray", upRegulated="red", downRegulated="green")
  gene <- rep("noDifference", nrow(volcanoSubstr))
  gene[logUp[logUp %in% withStat]] <- "upRegulated"
  gene[logDown[logDown %in% withStat]] <- "downRegulated"
  
  plot <- ggplot(data.frame(volcanoSubstr), aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(size=1.2) +
    geom_hline(yintercept=-log10(adjPValueThreshold), color="orange") +
    geom_vline(xintercept=-lfcThreshold, color="green") +
    geom_vline(xintercept=lfcThreshold, color="red") +
    aes(colour=gene) +
    scale_colour_manual(values=colours) +
    ggtitle("Volcano plot showing distribution of lfc vs adjp for expressed genes")
  
  print(plot)
}


candidateGeneExpressionBoxPlot <- function(geneOfInterestNormalExpression) {
  plot <- ggplot(geneOfInterestNormalExpression, aes(x=group, y=count, colour=group)) +
    geom_point(position=position_jitter(w=0.1,h=0), size=2) +
    scale_y_log10() +
    ggtitle("Distribution of normalised read counts for experimental conditions") +
    xlab("Experimental Condition") +
    ylab("normalised read count - log10 scale") +
    scale_colour_brewer(palette="Paired")
  print(plot)
}



pcaPlot <- function() {
  pcaMatrix <- counts(deSeqObj)
  md <- prep(t(pcaMatrix), scale="none", center=TRUE)
  pca <- pca(md, method="svd", center=TRUE, nPcs=3)
  xdata <- as.data.frame(pca@scores)
  xdata <- cbind(xdata, group=studyDesign$group[match(rownames(xdata), rownames(studyDesign))])
  x <- 1
  y <- 2
  xpercent <- round(pca@R2[x]*100, digits=1)
  ypercent <- round(pca@R2[y]*100, digits=1)
  xlab <- paste("Prin.Comp. ",x," (",xpercent,"%)",sep="")
  ylab <- paste("Prin.Comp. ",y," (",ypercent,"%)",sep="")
  
  plot <- ggplot(xdata, aes(x=PC1, y=PC2, colour=group)) +
    geom_point(size=5, shape=18) +
    scale_colour_brewer(palette="Paired") +
    geom_vline(xintercept=0, color="darkgray") +
    geom_hline(yintercept=0, color="darkgray") +
    ggtitle("PCA analysis of experimental samples") +
    ylab(paste("PC2 (",ypercent,"%)",sep=""))  +
    xlab(paste("PC1 (",xpercent,"%)",sep=""))
  
  print(plot)
}




