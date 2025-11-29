## Set working directory to the folder containing the count matrix and sample metadata
setwd("/path/to/project/diff_sites_fev23")
library("edgeR")
library("limma")
library("writexl")
library(pheatmap)
library(sva)
library(randomcoloR)
library(gtools)
library(org.Hs.eg.db)
library(Biobase)
library(genefilter)
library(tximport)
library(AnnotationDbi)
library(ggbiplot)

#txdb <- loadDb(file = "/path/to/bioinfo/R/txdb.gencode39.sqlite")
## Read gene-level count matrix (rows = genes, columns = samples) produced by the RNA-seq pipeline
counts <- read.table("counts_hisat_refSeq.tsv", header=1, sep="\t",row.names = 1)

## Drop auxiliary annotation column from the count matrix, keeping only sample columns
counts <- counts[,-c(1)]

## Read sample metadata (phenotype/clinical information) for each RNA-seq sample
pData <- read.table("pData.txt", header=1, sep="\t",row.names = 1)

## Ensure age is numeric so it can be used in exploratory analyses if needed
pData$Age <- as.numeric(pData$Age)

## Build design matrix with one coefficient per brain site; used for filtering and later models
des <- model.matrix(~0+Site, data=pData)
colnames(des)[1:length(levels(factor(pData$Site)))] <- levels(factor(pData$Site))

## Create DGEList object, normalise library sizes (TMM) and filter lowly expressed genes
## This improves variance estimates and reduces the multiple testing burden
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge, method = "TMM")
k <- filterByExpr(dge, design=des)
dge <- dge[k,]

#sy2 <- mapIds(org.Hs.eg.db, keys=rownames(dge), column="SYMBOL", keytype="ENTREZID")
sy2 <- rownames(dge)
#dge <- dge[!is.na(sy2),]

## Apply voom transformation to obtain log2-CPM values with precision weights for limma models
v <- voom(dge,plot=T)


## removing unwanted variation
## Estimate surrogate variables (SVs) capturing unmodelled technical variation (e.g. batch effects)
des0 <- model.matrix(~1, data=pData)
sva.res <- sva(v$E, des,des0)

## Extract SVs and append them to the sample metadata for later use in the design matrix
cova <- sva.res$sv
colnames(cova) <- paste("sva",1:ncol(cova),sep="")
pData <- cbind(pData,cova)
rownames(pData) <- colnames(v$E)

## Create batch-corrected expression matrix for visualisation, keeping the site effect in the design
matr <- removeBatchEffect(v, covariates=cova, design=des)
#matr <- v$E

## MDS plot: check if samples cluster mainly by brain site after normalisation and batch correction
pdf(file="MDSplot.pdf")
plotMDS(matr,labels=as.character(pData$Site),main="diff. sites")
dev.off()

# Faceting Attributes

library(ggplot2)
library(gridExtra)
library(ggfortify)
library(plotly)
library(ggbiplot)

# Faceting Attributes
#matr <- v$E
## PCA on the batch-corrected matrix to inspect global structure and separation between sites
pc.c <- prcomp(matr)

## PCA on samples (transpose) to use in biplots coloured by brain site
pc <- prcomp(t(matr))
pci <- data.frame(pc$x, Type=pData$Site)
percentage <- summary(pc)$importance[2,]*100
percentage <- paste( colnames(pci), "(", paste( as.character(percentage), "%", ")", sep="") )

#ko <- randomColor(length(levels(factor(pData$Site))), hue="random")

## Manual colour palette for brain sites to keep colours consistent across plots
ko <- c("#04416b", "#779ace", "#f9a7e1", "#ce6550", "#fff963", "#6aaa1b", "#73e2d7", "#5ef990")

## PCA biplots of different PC combinations to visualise how sites separate in low dimensions
p1 <- ggbiplot::ggbiplot(pc, scale.=T, groups=pData$Site, choices=c(1,2),
                         var.axes=F, ellipse=T, obs.scale = 1, var.scale = 1) +
  labs(title="diff. sites") +
  scale_color_manual(values=ko) +
  theme_bw()

p2 <- ggbiplot::ggbiplot(pc, scale.=T, groups=pData$Site, choices=c(1,3),
                         var.axes=F, ellipse=T, obs.scale = 1, var.scale = 1) +
  labs(title="diff. sites") +
  scale_color_manual(values=ko) +
  theme_bw() 

p3 <- ggbiplot::ggbiplot(pc, scale.=T, groups=pData$Site, choices=c(2,3),
                         var.axes=F, ellipse=T, obs.scale = 1, var.scale = 1) +
  labs(title="dm") +
  scale_color_manual(values=ko) +
  theme_bw() 

## Save PCA plots as PDFs for QC and potential figures
pdf(file="pca_plot_pc1_2.pdf")
p1
dev.off()

pdf(file="pca_plot2_pc1_3.pdf")
p2
dev.off()

pdf(file=".pca_plot2_pc2_3.pdf")
p3
dev.off()

## Final design matrix: brain sites of interest plus surrogate variables as covariates
## SVs adjust for unmeasured confounding while preserving the site coefficients
des.cor <- model.matrix(~0+Site+sva1+sva2+sva3+sva4+sva5+sva6+sva7+sva8+sva9+sva10+sva11, data=pData )
colnames(des.cor)[1:length(levels(factor(pData$Site)))] <- levels(factor(pData$Site))

## Fit limma linear model on voom-transformed data with the extended design
fit <- lmFit(v,des.cor)

## Build all pairwise contrasts between brain sites to test differences for every site combination
cmp <- combinations(n=length(levels(factor(pData$Site))), v=levels(factor(pData$Site)),r=2)
cm <- data.frame(matrix(nrow=length(colnames(des.cor)),ncol=nrow(cmp),data=0))
rownames(cm) <- colnames(des.cor)

for(i in 1:nrow(cmp)){
  cm[cmp[i,1],i] <- 1
  cm[cmp[i,2],i] <- -1
  colnames(cm)[i] <- paste(cmp[i,1],"-",cmp[i,2],sep="")
}
cm <- as.matrix(cm)
fit1 <- contrasts.fit(fit,cm)
fit1 <- eBayes(fit1) #, robust = T)

## Summarise significant genes under different log2FC and FDR cutoffs
## This allows us to examine DEGs under more or less stringent criteria
dt05 <- decideTests(fit1, adjust.method="BH", p.value=0.05, lfc=2)
dt1 <- decideTests(fit1, adjust.method="BH", p.value=0.05, lfc=1)
dt2 <- decideTests(fit1, adjust.method="BH", p.value=0.05, lfc=0)

dt05.1 <- decideTests(fit1, adjust.method="BH", p.value=0.1, lfc=2)
dt1.1 <- decideTests(fit1, adjust.method="BH", p.value=0.1, lfc=1)
dt2.1 <- decideTests(fit1, adjust.method="BH", p.value=0.1, lfc=0)

## Export DEG summaries (numbers of up/down genes per contrast) to an Excel workbook
library("openxlsx")
wb <-openxlsx::createWorkbook("dt")

addWorksheet(wb, "lfc2_p005")
writeData(wb, sheet=1, summary(dt05) )
addWorksheet(wb, "lfc1_p005")
writeData(wb, sheet=2, summary(dt1) )
addWorksheet(wb, "lfc0_p005")
writeData(wb, sheet=3, summary(dt2) )
addWorksheet(wb, "lfc2_p01")
writeData(wb, sheet=4, summary(dt05.1) )
addWorksheet(wb, "lfc1_p01")
writeData(wb, sheet=5, summary(dt1.1) )
addWorksheet(wb, "lfc0_p01")
writeData(wb, sheet=6, summary(dt2.1) )

openxlsx::saveWorkbook(wb,"DEG_results.xlsx",overwrite = T)

## Heatmaps of differentially expressed genes and downstream clustering
## Start with all genes significant at FDR < 0.05 (no log2FC cutoff) in at least one contrast

#1 - all

## Reorder dendrogram using the first singular vector to improve interpretability of heatmaps
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

## Select all genes with FDR < 0.05 across contrasts to build the expression matrix for clustering
res1.plot <- limma::topTable(fit1, number=Inf, p.value=0.05, adjust.method  ="BH", lfc=0)
mat2.plot <- matr[rownames(res1.plot),]

sy2 <- rownames(mat2.plot)

## Column annotation: brain site and sex for each sample
df <- data.frame(dm=pData$Site, sex=pData$Gender)
rownames(df) <- colnames(mat2.plot)

### SOM: self-organising map to group genes with similar expression patterns across samples
library("supraHex")
mat2.som <- mat2.plot
#rownames(mat2.som) <- unname(sy2)

## Centre each gene across samples before training the SOM to focus on relative patterns
matr_som <- mat2.som - matrix(rep(apply(mat2.som,1,mean), ncol(mat2.som)), ncol=ncol(mat2.som))
colnames(matr_som) <- pData$Site
#get trained

## recommended SOM size: https://bioconductor.org/packages/release/bioc/vignettes/oposSOM/inst/doc/Vignette.pdf
##                < 100 samples    100 - 500    500 - 1000     1000 - 5000     > 5000
## < 1000 genes      20 x 20        25 x 25      30 x 30         35 x 35       40 x 40
## 1000 - 10000      30 x 30        35 x 35      40 x 40         45 x 45       50 x 50
## 10000 - 50000     40 x 40        45 x 45      50 x 50         55 x 55       60 x 60

## Train SOM using supraHex; grid size chosen based on recommended settings for the number of genes/samples
sMap <- sPipeline(matr_som,xdim=40, ydim=40) #, algorithm="sequential",finetuneSustain = T )

## Visualise SOM component maps to inspect how expression patterns are distributed across the grid
visHexMulComp(sMap, title.rotate = 15, colormap = "jet", ncolors=20) ## click twice to open window
dev.list() ## choose quartz
dev.set(4) # 4 is the number of quartz
dev.copy(pdf, "SOM_map.pdf")
dev.off()

## Reorder SOM components based on correlation structure to aid interpretation
sReorder <- sCompReorder(matr_som,metric="pearson",amplifier=3, 
                         algorithm="sequential",neighKernel = "gaussian")

visCompReorder(sMap,sReorder,title.rotate=15,colormap = "jet", ncolors=20)
dev.list() ## choose quartz
dev.set(5) # 4 is the number of quartz
dev.copy(pdf, "SOM_map_reordered.pdf")
dev.off()

## Identify meta-clusters of SOM nodes and visualise them
sBase <- sDmatCluster(sMap=sMap) #,which_neigh = 1,
#distMeasure = "mean",clusterLinkage = "average",reindexSeed = "none")
#sBase <- sDmatCluster(sMap=sMap,which_neigh = 1,
#                      distMeasure = "max",clusterLinkage = "complete",reindexSeed = "none" )
visDmatCluster(sMap, sBase, colormap="jet")
dev.list() ## choose quartz
dev.set(6) # 4 is the number of quartz
dev.copy(pdf, "SOM_clusters.pdf")
dev.off()

## Assign each gene to a SOM cluster and export this mapping for downstream analyses
shuuu <- sWriteData(sMap,matr_som,sBase)
shuk <- shuuu[order(shuuu$Cluster_base, shuuu$Qerr_distance, decreasing=c(F,T)),]
write.csv(shuk, file="SOM_clusters.csv",quote = F, row.names=F)

shuk <- shuuu[order(shuuu$Cluster_base, shuuu$Qerr_distance, decreasing=c(F,T)),]
shuk <- shuk[!is.na(shuk$ID),]
rn <- shuk$ID
shuk <- data.frame(shuk)[,4,drop=F]
rownames(shuk) <- make.unique(rn)
colnames(shuk) <- "cluster"
shuk$cluster <- as.character(shuk$cluster)

mat3.plot <- matr[rownames(res1.plot),]
#rownames(mat3.plot) <- make.unique(mapIds(org.Rn.eg.db, keys=rownames(mat3.plot), column="SYMBOL", keytype="ENTREZID"))
mat3.plot <- mat3.plot[rownames(shuk),]

## Number of SOM clusters and colour palettes for cluster and site annotations
kk <- max(shuuu$Cluster_base)

kv <- randomColor(kk, hue="random")

ko.c1 <- ko
names(kv) <- as.character(1:kk)
names(ko.c1) <- levels(factor((pData$Site)))

ann_colors <- list(cluster=kv,
                   cond=ko.c1)


### reordering the matrix

## Order samples by brain site so heatmaps show clear site-specific blocks
pData.s <- pData[order(pData$Site),]

mat3.plot <- mat3.plot[,rownames(pData.s)]
df <- df[rownames(pData.s),,drop=F]

## Heatmap showing SOM-based gene clusters (rows) and samples ordered by brain site (columns)
pheatmap(mat3.plot, scale="row",  border_color = NA, ann_colors = ann_colors, annotation_row = shuk, 
         color=colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         annotation_col = df,
         show_colnames  = F, show_rownames = F,
         cluster_cols = F, cluster_rows = F, file = "SOM_heatmap_clusters_zoomout.pdf")

pheatmap(mat3.plot, scale="row",  border_color = NA, ann_colors = ann_colors, annotation_row = shuk, 
         color=colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         annotation_col = df,
         show_colnames  = F, show_rownames = T,h=1500,fontsize_row=6,
         cluster_cols = F, cluster_rows = F, file = "SOM_heatmap_clusters.pdf")


### GO and pathway enrichment for SOM clusters
### For each cluster, test enrichment in GO BP, KEGG and MSigDB gene sets
library("clusterProfiler")
library("openxlsx")
library("WebGestaltR")
library("msigdbr")
H <-  msigdbr(species = "human", category = "H")
C1 <- msigdbr(species = "human", category = "C1")
C2.CGP <- msigdbr(species = "human", category = "C2", subcategory = "CGP")
C2.CP <- msigdbr(species = "human", category = "C2", subcategory = "CP")
C3.TF <- msigdbr(species = "human", category = "C3", subcategory = "TFT:GTRD")


#rownames(mat2.plot) <- paste(names(sy2),";",sy2,sep="")
## Universe of genes: all genes present in the matrix, used as background for enrichment tests
univ <- rownames(matr)
write.table(univ, file="univ.txt",quote=F,col.names=F,row.names=F)
wb <-openxlsx::createWorkbook("enrichGO")
wc <-openxlsx::createWorkbook("enrichKEGG")
wa <-openxlsx::createWorkbook("msigdb")

all_sets <- msigdbr(species="Homo sapiens")
msigdbr_t2g = all_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

## Loop over clusters to run over-representation analysis and save results
for(i in 1:kk){
  enr <- enricher(shuuu$ID[shuuu$Cluster_base==i],pAdjustMethod="none", TERM2GENE=msigdbr_t2g)
  addWorksheet(wa, paste("cluster",i,sep=" "))
  writeData(wa, sheet=i, enr, colNames=T, rowNames=T)
}

openxlsx::saveWorkbook(wa,"enricher_MSIGDB_integrated_clusters.xlsx",overwrite = T)

for(i in 1:kk){
  clu <- enrichGO(shuuu$ID[shuuu$Cluster_base==i],pAdjustMethod="none",
                  'org.Hs.eg.db', ont="BP", qvalueCutoff=1, keyType = "SYMBOL",pvalueCutoff=0.1)
  addWorksheet(wb, paste("cluster",i,sep=" "))
  writeData(wb, sheet=i, clu, colNames=T, rowNames=T)
  if(dim(clu)[1]>0){
    pdf(paste("dotplot_cluster",i,"_integrated_GO.pdf",sep=""))
    try(print(dotplot(clu, orderBy = "x", title=paste("cluster",i,sep=" "),font.size=6)))
    dev.off()
  }
  
  #eee <- enrichr(shuuu$ID[shuuu$Cluster_base==i], dbs)
  
  entrezids <- mapIds(org.Hs.eg.db, keys=shuuu$ID[shuuu$Cluster_base==i], column="ENTREZID", keytype="SYMBOL")
  eee <- enrichKEGG(unname(entrezids), organism="hsa", keyType="kegg",pAdjustMethod = "none")
  addWorksheet(wc, paste("cluster",i,sep=" "))
  writeData(wc, sheet=i, eee, colNames=T, rowNames=T)
  if(!is.null(eee)){
    pdf(paste("dotplot_cluster",i,"_integrated_KEGG.pdf",sep=""))
    try(print(dotplot(eee, orderBy = "x", title=paste("cluster",i,sep=" "),font.size=6)))
    dev.off()
  }
  
  #if (!file.exists("./webgestalt")){
  #  dir.create("./webgestalt")
  #}
  #write.table(shuuu$ID[shuuu$Cluster_base==i],file="a.txt",quote=F,col.names=F,row.names = F )
  #if(!is.null(eee)){
  #  enres <- WebGestaltR(enrichMethod = "ORA", 
  #                     organism = "rnorvegicus",
  #                     sigMethod = "fdr",
  #                     fdrThr = 0.25,
  #                     enrichDatabase = "geneontology_Biological_Process_noRedundant",
  #                     interestGeneFile="a.txt",
  #                     interestGeneType="genesymbol",
  #                     referenceGeneFile="univ.txt", 
  #                     referenceGeneType = "genesymbol",
  #                     isOutput = T,
  #                     projectName = paste("GO_cluster_",i,"_integrated", sep=""),
  #                     outputDirectory = "./webgestalt")
  #  enres <- WebGestaltR(enrichMethod = "ORA", 
  #                     organism = "rnorvegicus",
  #                     sigMethod = "fdr",
  #                     fdrThr = 0.25,
  #                     enrichDatabase = "pathway_KEGG",
  #                     interestGeneFile="a.txt",
  #                     interestGeneType="genesymbol",
  #                     referenceGeneFile="univ.txt", 
  #                     referenceGeneType = "genesymbol",
  #                     isOutput = T,
  #                     projectName = paste("KEGG_cluster_",i,"_integrated", sep=""),
  #                     outputDirectory = "./webgestalt")
  #  enres <- WebGestaltR(enrichMethod = "ORA", 
  #                     organism = "rnorvegicus",
  #                     sigMethod = "fdr",
  #                     fdrThr = 0.25,
  #                     enrichDatabase = "network_Transcription_Factor_target",
  #                     interestGeneFile="a.txt",
  #                     interestGeneType="genesymbol",
  #                     referenceGeneFile="univ.txt", 
  #                     referenceGeneType = "genesymbol",
  #                     isOutput = T,
  #                     projectName = paste("TFT_cluster_",i,"_integrated", sep=""),
  #                     outputDirectory = "./webgestalt")
  
  # enres <- WebGestaltR(enrichMethod = "NTA", 
  #                     organism = "rnorvegicus",
  #                     sigMethod = "top",
  #                     topThr = 10,
  #                     enrichDatabase = "network_PPI_BIOGRID",
  #                     networkConstructionMethod = "Network_Expansion",
  #                     interestGeneFile="a.txt",
  #                     interestGeneType="genesymbol",
  #                     referenceGeneFile="univ.txt", 
  #                     referenceGeneType = "genesymbol",
  #                     isOutput = T,
  #                     NeighborNum = 10,
  #                     projectName = paste("PPI_BIOGRID_cluster_",i,"_integrated", sep=""),
  #                     outputDirectory = "./webgestalt")
  #}
} 

openxlsx::saveWorkbook(wb,"enrichGO_integrated_clusters.xlsx",overwrite = T)
openxlsx::saveWorkbook(wc,"enrichKEGG_integrated_clusters.xlsx",overwrite = T)

## DEG table: per-contrast statistics and mean expression by site for all genes

## Collect logFC, P.Value and adj.P.Val for each pairwise contrast into a single table
tabi <- data.frame(matrix(nrow=dim(matr)[1],ncol=0))
for(i in 1:length(colnames(fit1$contrasts))){
  tabi[,(dim(tabi)[2]+1):(dim(tabi)[2]+3)] <- topTable(fit1, coef=i, sort.by = "none", number=Inf, p.value=1, adjust.method  ="BH")[,c(1,4,5)]
  colnames(tabi)[(dim(tabi)[2]-2):(dim(tabi)[2])] <- c(paste(colnames(fit1$contrasts)[i],"_LogFC",sep=""),"P.Value","adj.P.Val")
}

## Compute average expression per site to help interpret the direction of change
meens <- data.frame(matrix(nrow=dim(matr)[1],ncol=0))
for(i in 1:length(colnames(des))){
  meens[,i] <- rowMeans(matr[,rownames(des)[which(des[,i,drop=F]==1)]])
  colnames(meens)[i] <- paste(colnames(des)[i],"_mean",sep="")
}

## Add gene-level annotation: symbol, Ensembl and Entrez IDs
symbol <- rownames(matr)
ensembl  <- mapIds(org.Hs.eg.db, keys=rownames(matr), column="ENSEMBL", keytype="SYMBOL")
entrez <-  mapIds(org.Hs.eg.db, keys=rownames(matr), column="ENTREZID", keytype="SYMBOL")

tab <- cbind(entrez,symbol,ensembl,matr,meens,tabi)
#tab2 <- cbind(ensembl,symbol,entrez,rpkm,meens,tabi)

## Write full gene-level expression and statistics table to Excel
wb <-openxlsx::createWorkbook("exp1")
addWorksheet(wb, "log_CPM")
writeData(wb, sheet=1, tab, colNames=T, rowNames=F)
#addWorksheet(wb, "rpkm")
#writeData(wb, sheet=2, tab2, colNames=T, rowNames=F)
openxlsx::saveWorkbook(wb,"diffsites_limma_sva.xlsx",overwrite = T)

library("GSVA")

library(GSEABase)

## Gene sets used for GSVA and ICA downstream analyses
## GMT files obtained from the MSigDB resource
## downloaded from https://www.gsea-msigdb.org/gsea/msigdb/mouse_geneset_resources.jsp

m1 <- getGmt("/path/to/GSVA_database/Homo/c1.all.v7.5.1.symbols.gmt.txt") ## positional
m2 <- getGmt("/path/to/GSVA_database/Homo/c2.all.v7.5.1.symbols.gmt.txt") ## curated
m3 <- getGmt("/path/to/GSVA_database/Homo/c3.all.v7.5.1.symbols.gmt.txt") ## regulatory
m4 <- getGmt("/path/to/GSVA_database/Homo/c4.all.v7.5.1.symbols.gmt.txt") ## Tumor ontoogy
m5 <- getGmt("/path/to/GSVA_database/Homo/c5.all.v7.5.1.symbols.gmt.txt") ## GO
m8 <- getGmt("/path/to/GSVA_database/Homo/c8.all.v7.5.1.symbols.gmt.txt") ## cell type
mh <- getGmt("/path/to/GSVA_database/Homo/h.all.v7.5.1.symbols.gmt.txt") ## hallmark

#sapply(strsplit(substr(rownames(mat2.plot),1,18)[res.clust$cluster==as.character(i)], split=";"), function(x) x[1])

## Build a custom gene set collection where each SOM cluster is treated as a gene set
cluj <- list()
for(i in 1:kk){
  cluj[[i]] <- shuuu$ID[shuuu$Cluster_base==i]
  names(cluj)[i] <- paste("cluster_",i,sep="")
}

matr_sy <- matr

GSVA_names <- c("m1","m2","m3","m4","m5","m8","mh","cluj")
list_names <- c("positional","curated","regulatory","TU_ontology",
                "GO","cell_type","hallmark","clusters")

cmp <- combinations(n=length(levels(factor(pData$Site))), v=levels(factor(pData$Site)),r=2)
cm <- data.frame(matrix(nrow=length(colnames(des)),ncol=nrow(cmp),data=0))
rownames(cm) <- colnames(des)

for(i in 1:nrow(cmp)){
  cm[cmp[i,1],i] <- 1
  cm[cmp[i,2],i] <- -1
  colnames(cm)[i] <- paste(cmp[i,1],"-",cmp[i,2],sep="")
}
cm <- as.matrix(cm)

## GSVA: compute gene set activity scores per sample for multiple MSigDB collections and SOM clusters
wb <-openxlsx::createWorkbook("GSVA")

## For each gene set collection, run GSVA, fit limma models and extract site contrasts at the gene set level
for(i in 1:length(GSVA_names)){
  re <- gsva(matr_sy, get(GSVA_names[i]), kcdf="Gaussian")
  ## Average GSVA scores per site to summarise gene set activity by brain region
  meenz <- data.frame(matrix(nrow=dim(re)[1],ncol=0))
  for(ji in 1:length(colnames(des))){
    meenz[,ji] <- rowMeans(re[,rownames(des)[which(des[,ji,drop=F]==1)]])
    colnames(meenz)[ji] <- paste(colnames(des)[ji],"_mean",sep="")
  }
  fit <- lmFit(re,des)
  fit1 <- contrasts.fit(fit,cm)
  fit1 <- eBayes(fit1) #, robust = T)
  tabi <- data.frame(matrix(nrow=dim(re)[1],ncol=0))
  for(j in 1:length(colnames(fit1$contrasts))){
    tabi[,(dim(tabi)[2]+1):(dim(tabi)[2]+3)] <- topTable(fit1, coef=j, sort.by = "none", number=Inf, p.value=1, adjust.method  ="BH")[,c(1,4,5)]
    colnames(tabi)[(dim(tabi)[2]-2):(dim(tabi)[2])] <- c(paste(colnames(fit1$contrasts)[j],"_LogFC",sep=""),"P.Value","adj.P.Val")
  }
  tab <- cbind(re,meenz,tabi)
  addWorksheet(wb, list_names[i])
  writeData(wb, sheet=i, tab, colNames=T, rowNames=T)
  ## For SOM-based cluster gene sets, also save boxplots of GSVA scores per site
  if(i==8){
    for(jj in 1:nrow(re)){
      pdf(paste("boxplot_GSVA_cluster_",jj,".pdf",sep=""))
      print(boxplot(re[jj,]~pData$Site, xlab="site",ylab="GSVA score", main=paste("GSVA_cluster_",jj,sep="")))
      dev.off()
    }
  }

    
}
         
openxlsx::saveWorkbook(wb,"GSVA_analysis.xlsx",overwrite = T)

## Cluster gene sets based on their GSVA profiles and visualise them in a heatmap
res <- pheatmap::pheatmap(re, annotation_col=df, show_colnames=F,show_rownames = F,clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",
                          scale="row",fontsize_row=6, clustering_callback = callback, clustering_method = "ward.D2",
                          cluster_cols =T, cluster_rows = T,color = colorRampPalette(c("navyblue", "white","firebrick3"))(50))
## Manually choose the number of gene set clusters (patterns) to cut the dendrogram
kk1=6
res.clust <- data.frame(cbind(re, pattern=cutree(res$tree_row, k=kk1)))
res.clust$pattern <- as.character(res.clust$pattern)

df_r <- data.frame(pattern=res.clust$pattern, stringsAsFactors = F)
rownames(df_r) <- rownames(re)


#write.csv(df_r, file = "./dm_noSVA/dm_clusters.csv")

## Colours for each GSVA pattern cluster and for sites
kp <- randomColor(kk1, hue="random")
#ko.c1 <- randomColor(length(levels(factor((pData$cond)))), hue="random")
ko.c1 <- ko
names(kp) <- as.character(1:kk1)
names(ko.c1) <- levels(factor((pData$Site)))

ann_colors <- list(pattern=kp,
                   cond=ko.c1)

pheatmap::pheatmap(re, annotation_col=df, annotation_row=df_r, show_colnames=F,clustering_distance_cols = "correlation",
                   show_rownames = T,clustering_distance_rows = "correlation", annotation_colors = ann_colors,file="GSVA_heatmap_clusters_zoomout.pdf",
                   scale="row",fontsize_row=6,main=paste("ward.D2, Pearson's correlation, k=",kk,sep=""), clustering_callback = callback, clustering_method = "ward.D2", treeheight_row = 150,
                   cluster_cols =T, cluster_rows = T,color = colorRampPalette(c("navyblue","white","firebrick3"))(50))

## add cluster info to orig matrix


### averaging  
for(i in 1:length(colnames(re))){
  re_avg <- data.frame(matrix(nrow=dim(re)[1],ncol=0))
  for(ji in 1:length(colnames(des))){
    re_avg[,ji] <- rowMeans(re[,rownames(des)[which(des[,ji,drop=F]==1)]])
    colnames(re_avg)[ji] <- paste(colnames(des)[ji],"_mean",sep="")
    rownames(re_avg) <- rownames(re)
  }
}

#re_avg <- re_avg[,c(2,4,5,3,1)]
#re_avg <- scale(re_avg, center = T, scale = F)
colnames(re_avg) <- c("A","B","C","D","E","F","G","H")

### kmeans clustering

## estimating the optimal number of clusters
#library(factoextra)
#fviz_nbclust(re_avg, FUN=kmeans, diss=dist(re_avg,"euclidean"), method = "wss") +
#  geom_vline(xintercept = 3, linetype = 2)

## optimal number: 4

#max_itr <- 10000
#n_clust  <-  4  ## number of cluster 
#set.seed(123) ## reproduce the cluster 
#kmeans_out  <- kmeans(re_avg,n_clust,iter.max = max_itr, nstart=100)

library(tidyverse)
library(ggrepel)
re_info <- data.frame(re_avg) %>% 
  mutate(clust = paste("pattern_", res.clust$pattern,sep = ""))
re_info$ori_cluster <- rownames(re_info)

names(kp) <- paste("pattern_",as.character(1:kk1),sep="")

## visualise  each cluster 
## Plot GSVA pattern profiles across sites for each cluster of gene sets
pp1 <- re_info %>% 
  gather(key = "variable" , value = "value", -c(9,10)) %>%  ### 1 is the index of column 'geneName' and 7 is the index of column 'clust'
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num)) + geom_point() +
  geom_line(alpha = 1, aes(col = as.character(clust))) + scale_color_manual(values=kp) +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
  facet_wrap(~clust) + scale_x_discrete(labels=c("A"="CER",
                                                 "B"="FRO",
                                                 "C"="HCP",
                                                 "D"="HTL",
                                                 "E"="OCC",
                                                 "F"="PAR",
                                                 "G"="SUB",
                                                 "H"="TEMP"))  


ggsave("GSVA_cluster_patterns.pdf", plot=pp1)


### ICA

library(ica)
library(msigdbr)

types <- c("H",paste("C",1:8,sep=""))

wd <- "./ICA"

ncomp <- 60
#lcpm <- matr 
## Run ICA (icaimax) on the expression matrix to identify independent transcriptional components
res <- icaimax(matr, ncomp, center=F, maxit=10000)
#res <- icajade(matr, ncomp, center=T, maxit=10000)
lcpm <- matr

i <- 1
j <- 1
glist <- list(NA)
ipi.l <- list(NA)
ini.l <- list(NA)

wk <-openxlsx::createWorkbook("gac")
jj=1
## For each component, identify genes with strong positive/negative loadings and save them
for(c in 1:ncomp){
  filp <- res$S[,c][res$S[,c] > 2.5]
  filn <- res$S[,c][res$S[,c] < -2.5]
  if(length(filp)>5){
    addWorksheet(wk, paste("Component_",c,"_positive_genes",sep=""))
    writeData(wk, sheet=jj, data.frame(symbol=names(filp)), colNames=T, rowNames=F)
    jj=jj+1
  }
  if(length(filn)>5){
  addWorksheet(wk, paste("Component_",c,"_negative_genes",sep=""))
  writeData(wk, sheet=jj, data.frame(symbol=names(filn)), colNames=T, rowNames=F)
  jj=jj+1
  }
}
openxlsx::saveWorkbook(wk,paste(wd,"/genes_and_components",ncomp,"_ICA.xlsx",sep=""),overwrite = T)

## For each ICA component, annotate leading positive/negative genes by MSigDB enrichment and GSVA
for(c in 1:ncomp){
  filp <- res$S[,c][res$S[,c] > 2.5]
  filn <- res$S[,c][res$S[,c] < -2.5]
  if(i==1){
    wb <-openxlsx::createWorkbook("ICA_msigdb")
    wc <-openxlsx::createWorkbook("ICA_gsva")
  }
  filp <- filp[order(filp, decreasing = T)]
  if(length(filp)>5){
    ipi <- gsva(lcpm,list(names(filp)[1:15]), Hs.H,min.sz=1, max.sz=Inf, verbose=T, kcdf="Gaussian")
    ipi.l[[c]] <- ipi
    names(ipi.l)[c] <-paste("Component_",c,"_positive",sep="") 
    addWorksheet(wc, paste("Component_",c,"_positive_GSVA",sep=""))
    writeData(wc, sheet=j, ipi, colNames=T, rowNames=T)
    j <- j+1
    #map <- matches(names(filp),ens.aggr$Group.1,all.x=T,all.y = F,list=F,indexes=T)
    #map <- map[order(map$x,decreasing=F),]
    #rnames <- ens.aggr$external_gene_name[map$y]
    #rnames[is.na(rnames)] <- names(filp)[is.na(map$y)]
    #rnames <- sapply(strsplit(rnames,";"),function(x){x[1]})
    syp <- names(filp)
    #filp <- filp[order(filp, decreasing = T)]
    #map <- match(names(filp),ens$ensembl_gene_id)
    #syp <- ens$external_gene_name[map]
    ## adding to RcisTarget object
    glist[[i]] <- syp
    names(glist)[i] <-paste("Component_",c,"_positive",sep="") 
    for(m in 1:length(types)){
      m_df <- msigdbr(species = "Homo sapiens", category=types[m]) %>% dplyr::select(gs_name, gene_symbol)
      alll <- enricher(syp, TERM2GENE=m_df)
      addWorksheet(wb, paste("Component_",c,"_positive_",types[m],sep=""))
      writeData(wb, sheet=i, alll, colNames=T, rowNames=T)
      i <- i+1
    }
    ## Heatmap of top positive genes for this component, ordered by GSVA score across samples
    mattp <- lcpm[names(filp),order(ipi)]
    rownames(mattp) <- syp
    fnn <- paste(wd,"/comp_",c,"_pos.pdf",sep="")
    pheatmap::pheatmap(mattp,scale="row",cluster_rows = F, cluster_cols = F,height=20,
                       annotation_names_col = F, filename = fnn, fontsize = 4,clustering_distance_rows = "correlation",
                       clustering_distance_cols = "correlation", annotation_col = df, annotation_colors = ann_colors,
                       color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  }
  filn <- filn[order(filn, decreasing = F)]
  if(length(filn)>5){
    ini <- gsva(lcpm,list(names(filn)[1:15]), Hs.H,min.sz=1, max.sz=Inf, verbose=T, kcdf="Gaussian")
    ini.l[[c]] <- ini
    names(ini.l)[c] <-paste("Component_",c,"_negative",sep="") 
    addWorksheet(wc, paste("Component_",c,"_negative_GSVA",sep=""))
    writeData(wc, sheet=j, ini, colNames=T, rowNames=T)
    j <- j+1
    #man <- matches(names(filn),ens.aggr$Group.1,all.x=T,all.y = F,list=F,indexes=T)
    #man <- man[order(man$x,decreasing=F),]
    #rnames <- ens.aggr$external_gene_name[man$y]
    #rnames[is.na(rnames)] <- names(filn)[is.na(man$y)]
    #rnames <- sapply(strsplit(rnames,";"),function(x){x[1]})
    syn <- names(filn)
    #filn <- filn[order(filn, decreasing = F)]
    #man <- match(names(filn),ens$ensembl_gene_id)
    #syn <- ens$external_gene_name[man]
    ## adding to RcisTarget object
    glist[[i]] <- syn
    names(glist)[i] <-paste("Component_",c,"_negative",sep="") 
    for(m in 1:length(types)){
      m_df <- msigdbr(species = "Homo sapiens", category=types[m]) %>% dplyr::select(gs_name, gene_symbol)
      alll <- enricher(syn, TERM2GENE=m_df)
      addWorksheet(wb, paste("Component_",c,"_negative_",types[m],sep=""))
      writeData(wb, sheet=i, alll, colNames=T, rowNames=T)
      i <- i+1
    }
    ## Heatmap of top negative genes for this component to visualise opposite expression patterns
    mattn <- lcpm[names(filn),order(ini)]
    rownames(mattn) <- syn
    fnn <- paste(wd,"/comp_",c,"._neg.pdf",sep="")
    pheatmap::pheatmap(mattn,scale="row",cluster_rows = F, cluster_cols = F,height=20,
                       annotation_names_col = F, filename = fnn, fontsize = 4,clustering_distance_rows = "correlation",
                       clustering_distance_cols = "correlation", annotation_col = df, annotation_colors = ann_colors,
                       color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  }
  
}

## Save all ICA-based enrichment and GSVA results for later inspection
openxlsx::saveWorkbook(wb,paste(wd,"/enricher_ICA_",ncomp,"_msigdb_ICA.xlsx",sep=""),overwrite = T)
openxlsx::saveWorkbook(wc,paste(wd,"/GSVA_ICA_",ncomp,".xlsx",sep=""),overwrite = T)
