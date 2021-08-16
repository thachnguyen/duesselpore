####################################################################################################
## Limma edgeR ######
####################

library(limma)
library(edgeR)

y <- DGEList(geneCounts, group=studyDesign$group, genes=rownames(geneCounts))
#options(digits=3)
y$samples

group <- factor(studyDesign$group)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
AveLogCPM <- aveLogCPM(y)

y <- calcNormFactors(y)

logCPM <- cpm(y, prior.count=2, log=TRUE)

colnames(logCPM) <- rownames(y$samples)

#Estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)

plotQLDisp(fit)
summary(fit$df.prior)

if (config$referenceGroup %in% levels(group)&(config$studyGroup %in% levels(group))){
  i <- match(config$referenceGroup, levels(group))
  j <-match(config$studyGroup, levels(group))
}else{
  i <-1
j <-2}

group1vs2<-makeContrasts(paste(levels(group)[j], levels(group)[i], sep = '-'), levels=design)

res <- glmQLFTest(fit, contrast=group1vs2)
topTags(res)
tr <- glmTreat(fit, contrast=group1vs2, lfc=log2(1.5))
topTags(tr)

o <- order(tr$table$PValue)

logCPM <- logCPM - rowMeans(logCPM)
logCPM <- logCPM[o[1:config$NumberOfTopGene],]

anno <- as.data.frame(studyDesign[, c("group","replicate")])
anno$replicate <- factor(anno$replicate)

variance_heatmap <- pheatmap(logCPM, 
                             annotation_col = anno, 
                             cluster_cols=config$cluster_col,
                             clustering_distance_cols = "correlation",
                             labels_row = mapIds(refdb, keys = substr(rownames(logCPM),1,15), 
                                                 column = "SYMBOL", keytype = "GENEID", multiVals = "first"),
                             labels_col = c(rownames(studyDesign)), drop_levels = TRUE, filename='Analysis/Results/heatmap.pdf')


organism <- org.Hs.eg.db
res_group01_group02 <- tr

res_group01_group02.filtered <- res_group01_group02 %>%
  as.data.frame() %>%
  dplyr::filter(abs(logFC) > config$lfcThreshold)

geneList <- res_group01_group02.filtered$logFC
names(geneList) <- res_group01_group02.filtered$genes
geneList <- sort(geneList, decreasing = TRUE)

options(repr.plot.width=12, repr.plot.height=15)

gse <- gseGO(geneList = geneList,
             ont = "ALL",
             keyType = "ENSEMBL",
             nPerm = 10000,
             #minGSSize = 5,
             #maxGSSize = 500,
             #pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "none",
             eps =0.0,
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
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
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
pathway_file <- list.files('.','.png$')
file.copy(pathway_file, 'Analysis/Results/')