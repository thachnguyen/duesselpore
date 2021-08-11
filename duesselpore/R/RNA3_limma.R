#########################################################
# Analyzing the differential expression using limma
# Thresholds have been defined in the config.yaml-File

library(limma)
library(edgeR)

y <- DGEList(geneCounts, group=studyDesign$group, genes=rownames(geneCounts))
options(digits=3)
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

colnames(logCPM) <- paste(y$samples$samples)

#Estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)

logCPM <- logCPM[o[1:config$NumberOfTopGene],]

plotQLDisp(fit)
summary(fit$df.prior)
group1vs2<-makeContrasts(paste(levels(group)[2], levels(group)[1], sep = '-'), levels=design)

res <- glmQLFTest(fit, contrast=group1vs2)
topTags(res)
tr <- glmTreat(fit, contrast=group1vs2, lfc=log2(1.5))
topTags(tr)

logCPM <- cpm(y, prior.count=2, log=TRUE)
colnames(logCPM) <- paste(y$samples$group, 1:3, sep="-")
o <- order(tr$table$PValue)

logCPM <- logCPM - rowMeans(logCPM)
logCPM <- logCPM[o[1:config$NumberOfTopGene],]

variance_heatmap <- pheatmap(logCPM, 
               annotation_col = anno, 
               cluster_cols=config$cluster_col,
               labels_row = mapIds(EnsDb.Hsapiens.v86, keys = substr(rownames(logCPM),1,15), 
                                   column = "SYMBOL", keytype = "GENEID", multiVals = "first"),
              labels_col = c(rownames(studyDesign)), drop_levels = FALSE, filename='Analysis/Results/heatmap.pdf')



organism <- org.Hs.eg.db
res_group01_group02 <- tr

res_group01_group02.filtered <- res_group01_group02 %>%
  as.data.frame() %>%
  dplyr::filter(abs(logFC) > 1.5)

geneList <- res_group01_group02.filtered$logFC
names(geneList) <- res_group01_group02.filtered$genes
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

