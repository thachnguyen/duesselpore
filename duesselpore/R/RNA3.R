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
         clustering_distance_cols = "correlation",
         labels_row = mapIds(refdb, keys = substr(rownames(mat),1,15), 
                             column = "SYMBOL", keytype = "GENEID", multiVals = "first"),
         labels_col = c(rownames(studyDesign)), drop_levels = TRUE, filename='Analysis/Results/heatmap.pdf')

organism <- org.Hs.eg.db

#####################
group <-factor(studyDesign$group)
if (config$referenceGroup %in% studyDesign$group&(config$studyGroup %in% studyDesign$group)){
  res_group01_group02 <- results(dds, contrast = c("group", config$studyGroup, config$referenceGroup))
}else{
  res_group01_group02 <- results(dds, contrast = c("group", levels(group)[2], levels(group)[1]))}
#####################


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

geneList1<-geneList
de_ens <- geneList1[abs(geneList1) > 1]
de <-as.data.frame(de_ens)

##################
#ensembl to entrez
entr1 = mapIds(org.Hs.eg.db,
               keys=row.names(de), 
               column="ENTREZID",
               keytype="ENSEMBL",
               multiVals="first")
entr1 <-as.data.frame(entr1)
de$entrez <-entr1$entr1
de <-na.omit(de)

de1 <- de$de_ens
names(de1) <-de$entrez
#################

edo <- enrichDGN(names(de1))
barplot(edo, showCategory=20)
ggsave("Analysis/Results/bar_plot_enrichDGN.pdf", width = 10, height = 8)

############
options(repr.plot.width=40, repr.plot.height=12)
edox <- setReadable(edo, 'org.Hs.eg.db', keyType = 'ENTREZID')
p1 <- cnetplot(edox, foldChange=de1)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=de1)
p3 <- cnetplot(edox, foldChange=de1, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p2, p3, ncol=2, labels=LETTERS[1:3], rel_widths=c(.8, 1.2))
#cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

cowplot::ggsave2("Analysis/Results/cnet_plot1.pdf", width = 20, height = 8)

p1 <- cnetplot(edox, node_label="category") 
p2 <- cnetplot(edox, node_label="gene") 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
cowplot::ggsave2("Analysis/Results/cnet_plot2.pdf", width = 20, height = 14)

p1 <- heatplot(edox)
p2 <- heatplot(edox, foldChange=de1)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
cowplot::ggsave2("Analysis/Results/heatmap_functional_classification.pdf", width = 20, height = 14)

hsa05321 <- pathview(gene.data  = de1,
                     pathway.id = config$pathway_ID,
                     species    = "hsa",
                     min.nnodes=3,
                     expand.node=TRUE,
                     )
pathway_file <- list.files('.','.png$')[1]
file.copy(pathway_file, 'Analysis/Results/pathway.png')


}
