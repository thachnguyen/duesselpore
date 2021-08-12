#########################################################
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
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), n = config$NumberOfTopGene)
mat <- assay(vsd)[topVarGenes,]
mat <- mat - rowMeans(mat)

mat_breaks <- seq(min(mat), max(mat), length.out = 10)
anno <- as.data.frame(colData(vsd)[, c("group","replicate")])
anno$replicate <- factor(anno$replicate) 
variance_heatmap <- pheatmap(mat, 
         annotation_col = anno, 
         cluster_cols=config$cluster_col,
         labels_row = mapIds(refdb, keys = substr(rownames(mat),1,15), 
                             column = "SYMBOL", keytype = "GENEID", multiVals = "first"),
         labels_col = c(rownames(studyDesign)), drop_levels = TRUE, filename='Analysis/Results/heatmap.pdf')

organism <- org.Hs.eg.db
res_group01_group02 <- results(dds)

res_group01_group02.filtered <- res_group01_group02 %>%
  as.data.frame() %>%
  dplyr::filter(abs(log2FoldChange) > config$lfcThreshold)

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

############
options(repr.plot.width=40, repr.plot.height=12)
edox <- setReadable(edo, 'org.Hs.eg.db', keyType = 'auto')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

cowplot::ggsave2("Analysis/Results/cnet_plot1.pdf", width = 20, height = 10)

p1 <- cnetplot(edox, node_label="category") 
p2 <- cnetplot(edox, node_label="gene") 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
cowplot::ggsave2("Analysis/Results/cnet_plot2.pdf", width = 20, height = 14)

p1 <- heatplot(edox)
p2 <- heatplot(edox, foldChange=geneList)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
cowplot::ggsave2("Analysis/Results/heatmap_functional_classification.pdf", width = 20, height = 14)

library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     kegg.dir="Analysis/Results/",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

}
