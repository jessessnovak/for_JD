library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(dplyr)
library(tidyverse)
library(Seurat)
library(patchwork)
#BiocManager::install("ComplexHeatmap") 
library(ComplexHeatmap)
library(RColorBrewer)
library(Cairo)
#install.packages("Cairo")
library(reticulate)
library(ggplot2)
#install.packages("leiden")
library("leiden")
library("igraph")
library("biomaRt")
#install.packages("EnhancedVolcano")
#library(EnhancedVolcano)

##################### experiment explanation ##################### 

#   3 time points: unwounded (UW), 3 days post wounding (D3),
#   and 5 days post wounding (D5)
#   2 genotypes: WT, and KO for each time point for a total of 6 samples
#   CMO301 = UW_WT, CMO302 = UW_KO
#   CMO303 = D3_WT, CMO304 = D3_KO
#   CMO305 = D5_WT, CMO306 = D5_KO
#   KOs heal wounds faster than WT
#   KOs may differentiate more than WT

##################### init #####################



#setwd("~/Dropbox (Personal)/Weill Cornell/PhD Years/Fuchs Lab/Fuchs Lab Notebook/Bioinformatics/FuchsLab/Phgdh_cKO_scRNAseq")
setwd("~/Dropbox (Personal)/Weill Cornell/PhD Years/Fuchs Lab/Fuchs Lab Notebook/Bioinformatics/FuchsLab/Phgdh_cKO_scRNAseq/PhgdhfloxedDremelUW_multiCMO2")
#setwd("~/txg-macos-v1.3.0/PhgdhfloxedDremelUW_multiCMO2")
data.dir = "raw_feature_bc_matrix"

#read in multiplexing metadata
multi.data <- read.csv("multiplexing_analysis/assignment_confidence_table.csv")
#make barcodes rownames 
multi.data <- relocate(multi.data, Assignment) %>% 
  relocate(Barcode)
rownames(multi.data) <- multi.data[,1]

#read in 10x counts and HTO data
hfsc.data <- Read10X(data.dir)
hfsc.data

# Load in the UMI matrix
hfsc.umis <- hfsc.data$`Gene Expression`

# Setup Seurat object with barcode assignments
hfsc.hashtag <- CreateSeuratObject(counts = hfsc.umis, 
                                   meta.data = multi.data)
#hfsc.hashtag@meta.data

#remove blank barcodes
hfsc.hashtag = subset(x = hfsc.hashtag, 
                      subset = Assignment != 'Blank')
#remove multiplets
hfsc.hashtag = subset(x = hfsc.hashtag, 
                      subset = Assignment != 'Multiplet')
#remove unassigned cells
hfsc.hashtag = subset(x = hfsc.hashtag, 
                      subset = Assignment != 'Unassigned')

hfsc.hashtag[["percent.mt"]] <- PercentageFeatureSet(hfsc.hashtag, pattern = "^mt-") #add percentage mito reads in metadata

#hfsc.hashtag.unfiltered <- hfsc.hashtag
#hfsc.hashtag <- hfsc.hashtag.unfiltered #revert


VlnPlot(hfsc.hashtag, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) #for some reason requires Cairo, which requires XQuartz

plot1 <- FeatureScatter(hfsc.hashtag, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hfsc.hashtag, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot1
plot2

summary(hfsc.hashtag@meta.data$nCount_RNA)
summary(hfsc.hashtag@meta.data$percent.mt)

hfsc.hashtag <- subset(hfsc.hashtag, subset = nFeature_RNA > 200 & percent.mt < 10)

hfsc.hashtag <- NormalizeData(hfsc.hashtag, normalization.method = "LogNormalize", scale.factor = 10000)
hfsc.hashtag <- FindVariableFeatures(hfsc.hashtag, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(hfsc.hashtag), 10) #not able to identify any

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(hfsc.hashtag)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#pre-processing to dim reduction
all.genes <- rownames(hfsc.hashtag)
hfsc.hashtag <- ScaleData(hfsc.hashtag, features = all.genes)
#pca dim reduction
hfsc.hashtag <- RunPCA(hfsc.hashtag, features = VariableFeatures(object = hfsc.hashtag))

#print(hfsc.hashtag[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(hfsc.hashtag, dims = 1:2, reduction = "pca")
##DimPlot(hfsc.hashtag, reduction = "pca")

#DimHeatmap(hfsc.hashtag, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(hfsc.hashtag, dims = 1:15, cells = 500, balanced = TRUE)

# use JackStraw to determine # PCs
#hfsc.hashtag <- JackStraw(hfsc.hashtag, num.replicate = 100)
#hfsc.hashtag <- ScoreJackStraw(hfsc.hashtag, dims = 1:20)
#JackStrawPlot(hfsc.hashtag, dims = 1:20)

ElbowPlot(hfsc.hashtag, 50) # determine # PCs

####################################Cluster Analysis####################################
hfsc.hashtag <- FindNeighbors(hfsc.hashtag)
hfsc.hashtag <- FindClusters(hfsc.hashtag, resolution = 0.8, 
                             algorithm = 4) #use leiden clustering
hfsc.hashtag <- RunUMAP(hfsc.hashtag, dims = 1:22)
DimPlot(hfsc.hashtag, reduction = "umap", label = T, 
        label.size = 4, label.box = T)

#FeaturePlot(hfsc.hashtag, features = "Areg", pt.size = 0.5)

#filter out most contaminating non-epithelial cells
hfsc.hashtag.filtered <- subset(x = hfsc.hashtag, subset = Ptprc < 0.05)
hfsc.hashtag.filtered <- subset(x = hfsc.hashtag.filtered, subset = Kit < 0.05)
hfsc.hashtag.filtered <- subset(x = hfsc.hashtag.filtered, subset = Pdgfra < 0.05)
hfsc.hashtag.filtered <- subset(x = hfsc.hashtag.filtered, subset = Pecam1 < 0.05)
#FeaturePlot(hfsc.hashtag.filtered, features = "Pecam1", pt.size = 0.5)
hfsc.hashtag <- hfsc.hashtag.filtered

#pick which identity
Idents(hfsc.hashtag) <- "seurat_clusters" #show timepoints for cells
#Idents(hfsc.hashtag) <- "Assignment" #show barcodes for cells
#Idents(mydata) <- "Phase" # show cell cycle state, not yet performed
#Idents(hfsc.hashtag) <- factor(x = Idents(hfsc.hashtag), levels = sort(levels(hfsc.hashtag)))
#DimPlot(hfsc.hashtag, reduction = "umap")

#FeaturePlot(hfsc.hashtag, features = c("Krt24", "Phgdh"), pt.size = 0.5)

#find positive markers for each cluster
hfsc.hashtag <- SetIdent(hfsc.hashtag, value = "seurat_clusters")
clust.markers <- FindAllMarkers(hfsc.hashtag, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(clust.markers, file = "all_cluster_markers.csv")

#top 20 per cluster
topG <- clust.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
head(topG)
DoHeatmap(hfsc.hashtag, features = topG$gene) + NoLegend()

write.csv(topG, file = "topgenes_clusters.csv")

#rename clusters with working names
new.cluster.ids <- c("HFSC_1", "Proliferative_1", "Proliferative_2", 
                     "Differentiating_1", "Wounded_Unwounded_HF", "Stressed_1",
                     "Lineage_Infidelity_1", "Migratory", "HFSC_2", 
                     "Hypoxic", "Stressed_2", "Upper_suprabasal",
                     "Unknown", "Cornifying")
names(new.cluster.ids) <- levels(hfsc.hashtag)
hfsc.hashtag <- RenameIdents(hfsc.hashtag, new.cluster.ids)
Idents(hfsc.hashtag) <- factor(x = Idents(hfsc.hashtag), levels = sort(levels(hfsc.hashtag)))
DimPlot(hfsc.hashtag, reduction = "umap", label = T, repel = T, 
        pt.size = 1, label.size = 5, label.box = T, label.color = "white") #+ NoLegend()


#saveRDS(hfsc.hashtag, "hfsc_hashtag_seurat.rds")

hfsc.hashtag <- readRDS("hfsc_hashtag_seurat.rds")

#VlnPlot(hfsc.hashtag, features = c("Nqo1", "Gclc","Epgn", 
                             "Hmox1", "Gpx2"))

#FeaturePlot(hfsc.hashtag, features = c("Sox9", "Klf5"), pt.size = 1.5)

################ recluster UW only ################ 
Idents(hfsc.hashtag) <- "Assignment" 
cell_values <- c("CMO301", "CMO302")
uw.hfsc <- subset(hfsc.hashtag, idents = cell_values, invert = FALSE)
uw.hfsc <- NormalizeData(uw.hfsc, normalization.method = "LogNormalize", scale.factor = 10000)
uw.hfsc <- FindVariableFeatures(uw.hfsc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(uw.hfsc), 10) 

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(uw.hfsc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#pre-processing to dim reduction
all.genes <- rownames(uw.hfsc)
uw.hfsc <- ScaleData(uw.hfsc, features = all.genes)
#pca dim reduction
uw.hfsc <- RunPCA(uw.hfsc, features = VariableFeatures(object = uw.hfsc))

#print(uw.hfsc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(uw.hfsc, dims = 1:2, reduction = "pca")
##DimPlot(uw.hfsc, reduction = "pca")

#DimHeatmap(uw.hfsc, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(uw.hfsc, dims = 1:15, cells = 500, balanced = TRUE)

# use JackStraw to determine # PCs
#uw.hfsc <- JackStraw(uw.hfsc, num.replicate = 100)
#uw.hfsc <- ScoreJackStraw(uw.hfsc, dims = 1:20)
#JackStrawPlot(uw.hfsc, dims = 1:20)

ElbowPlot(uw.hfsc, 50) # determine # PCs

####################################Cluster Analysis####################################
uw.hfsc <- FindNeighbors(uw.hfsc)
uw.hfsc <- FindClusters(uw.hfsc, resolution = 0.8, algorithm = 4) #use leiden clustering
#head(Idents(uw.hfsc), 5)
uw.hfsc <- RunUMAP(uw.hfsc, dims = 1:16)
DimPlot(uw.hfsc, reduction = "umap", label = T, label.size = 4, label.box = T)

Idents(uw.hfsc) <- "seurat_clusters" #show timepoints for cells
Idents(uw.hfsc) <- "Assignment" #show barcodes for cells
#Idents(uw.hfsc) <- "Phase" # show cell cycle state
Idents(uw.hfsc) <- factor(x = Idents(uw.hfsc), levels = sort(levels(uw.hfsc)))
DimPlot(uw.hfsc, reduction = "umap", label = T, label.size = 4, label.box = T)

Idents(uw.hfsc) <- "Assignment" #show barcodes for cells
Idents(hfsc.hashtag) <- factor(x = Idents(hfsc.hashtag), levels = sort(levels(hfsc.hashtag)))
HFSC_uw_Bimod <- FindMarkers(uw.hfsc, ident.1 = "CMO301", ident.2 = "CMO302", logfc.threshold = 0.25, test.use = "bimod")
head(HFSC_uw_Bimod)
write.csv(HFSC_uw_Bimod, "HFSC_uw_Bimod.csv")

VlnPlot(uw.hfsc, features = "Lrch1", pt.size = 0.5)

FeaturePlot(uw.hfsc, features = "Phgdh", pt.size = 0.5)
################ recluster D3 only ################ 
Idents(hfsc.hashtag) <- "Assignment" 
Idents(uw.hfsc) <- factor(x = Idents(uw.hfsc), levels = sort(levels(uw.hfsc)))
cell_values <- c("CMO303", "CMO304")
d3.hfsc <- subset(hfsc.hashtag, idents = cell_values, invert = FALSE)
d3.hfsc
d3.hfsc <- NormalizeData(d3.hfsc, normalization.method = "LogNormalize", scale.factor = 10000)
d3.hfsc <- FindVariableFeatures(d3.hfsc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(d3.hfsc), 10) 

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(d3.hfsc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#pre-processing to dim reduction
all.genes <- rownames(d3.hfsc)
d3.hfsc <- ScaleData(d3.hfsc, features = all.genes)
#pca dim reduction
d3.hfsc <- RunPCA(d3.hfsc, features = VariableFeatures(object = d3.hfsc))

#print(d3.hfsc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(d3.hfsc, dims = 1:2, reduction = "pca")
##DimPlot(d3.hfsc, reduction = "pca")

#DimHeatmap(d3.hfsc, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(d3.hfsc, dims = 1:15, cells = 500, balanced = TRUE)

# use JackStraw to determine # PCs
#d3.hfsc <- JackStraw(d3.hfsc, num.replicate = 100)
#d3.hfsc <- ScoreJackStraw(d3.hfsc, dims = 1:20)
#JackStrawPlot(d3.hfsc, dims = 1:20)

ElbowPlot(d3.hfsc, 50) # determine # PCs

####################################Cluster Analysis####################################
d3.hfsc <- FindNeighbors(d3.hfsc)
d3.hfsc <- FindClusters(d3.hfsc, resolution = 0.8, algorithm = 4) #use leiden clustering
#head(Idents(d3.hfsc), 5)
d3.hfsc <- RunUMAP(d3.hfsc, dims = 1:19)
DimPlot(d3.hfsc, reduction = "umap", label = T, label.size = 4, label.box = T)

Idents(d3.hfsc) <- "seurat_clusters" #show timepoints for cells
Idents(d3.hfsc) <- "Assignment" #show barcodes for cells
#Idents(d3.hfsc) <- "Phase" # show cell cycle state
Idents(d3.hfsc) <- factor(x = Idents(d3.hfsc), levels = sort(levels(d3.hfsc)))
DimPlot(d3.hfsc, reduction = "umap", label = T, label.size = 4, label.box = T)

HFSC_d3_Bimod <- FindMarkers(d3.hfsc, ident.1 = "CMO303", ident.2 = "CMO304", logfc.threshold = 0.25, test.use = "bimod")
head(HFSC_d3_Bimod)
write.csv(HFSC_d3_Bimod, "HFSC_d3_Bimod.csv")
VlnPlot(d3.hfsc, features = "Phgdh", pt.size = 0.5)

FeaturePlot(d3.hfsc, features = "Cdh3", pt.size = 0.5)

################ recluster D5 only ################ 
Idents(hfsc.hashtag) <- "Assignment" 
Idents(hfsc.hashtag) <- factor(x = Idents(hfsc.hashtag), levels = sort(levels(hfsc.hashtag)))
cell_values <- c("CMO305", "CMO306")
d5.hfsc <- subset(hfsc.hashtag, idents = cell_values, invert = FALSE)
d5.hfsc
d5.hfsc <- NormalizeData(d5.hfsc, normalization.method = "LogNormalize", scale.factor = 10000)
d5.hfsc <- FindVariableFeatures(d5.hfsc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(d5.hfsc), 10) 

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(d5.hfsc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#pre-processing to dim reduction
all.genes <- rownames(d5.hfsc)
d5.hfsc <- ScaleData(d5.hfsc, features = all.genes)
#pca dim reduction
d5.hfsc <- RunPCA(d5.hfsc, features = VariableFeatures(object = d5.hfsc))

#print(d5.hfsc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(d5.hfsc, dims = 1:2, reduction = "pca")
##DimPlot(d5.hfsc, reduction = "pca")

#DimHeatmap(d5.hfsc, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(d5.hfsc, dims = 1:15, cells = 500, balanced = TRUE)

# use JackStraw to determine # PCs
#d5.hfsc <- JackStraw(d5.hfsc, num.replicate = 100)
#d5.hfsc <- ScoreJackStraw(d5.hfsc, dims = 1:20)
#JackStrawPlot(d5.hfsc, dims = 1:20)

ElbowPlot(d5.hfsc, 50) # determine # PCs

####################################Cluster Analysis####################################
d5.hfsc <- FindNeighbors(d5.hfsc)
d5.hfsc <- FindClusters(d5.hfsc, resolution = 0.8, algorithm = 4) #use leiden clustering
#head(Idents(d5.hfsc), 5)
d5.hfsc <- RunUMAP(d5.hfsc, dims = 1:22)
DimPlot(d5.hfsc, reduction = "umap", label = T, label.size = 4, label.box = T)

VlnPlot(d5.hfsc, features = "Sox5", pt.size = 0.5)
?VlnPlot()

#Idents(d5.hfsc) <- "seurat_clusters" #show timepoints for cells
Idents(d5.hfsc) <- "Assignment" #show barcodes for cells
#Idents(d5.hfsc) <- "Phase" # show cell cycle state
Idents(d5.hfsc) <- factor(x = Idents(d5.hfsc), levels = sort(levels(d5.hfsc)))
DimPlot(d5.hfsc, reduction = "umap", label = T, label.size = 4, label.box = T)

HFSC_d5_Bimod <- FindMarkers(d5.hfsc, ident.1 = "CMO305", ident.2 = "CMO306", logfc.threshold = 0.25, test.use = "bimod")
HFSC_d5_Bimod
write.csv(HFSC_d5_Bimod, "HFSC_d5_Bimod.csv")
VlnPlot(d5.hfsc, features = "Ass1", pt.size = 0.5)

FeaturePlot(d5.hfsc, features = "Itgb1", pt.size = 0.5)

##################### SLINGSHOT ##################### 

#topology
remotes::install_github("kstreet13/bioc2020trajectories")
sce <- as.SingleCellExperiment(hfsc.hashtag, assay = 'RNA')
shuffle <- sample(ncol(sce))
layout(matrix(1:2, nrow = 1))
par(mar = c(4.5,4,1,1))

plot(reducedDims(sce)$UMAP[shuffle, ],
     asp = 4, pch = 16, xlab = 'UMAP-1', ylab = 'UMAP-2',
     col = alpha(c(1:2)[factor(colData(sce)$Assignment)][shuffle], alpha = 0.5))


scores <- bioc2020trajectories::imbalance_score(rd = reducedDims(sce)$UMAP,
  cl = colData(sce)$Assignment,
  k = 20, smooth = 40)

grad <- viridis::viridis(10, begin = 0, end = 1)
names(grad) <- levels(cut(scores$scaled_scores, breaks = 10))
plot(reducedDims(sce)$UMAP, col = grad[cut(scores$scaled_scores, breaks = 10)],
     asp = 2, pch = 16, xlab = 'UMAP-1', ylab = 'UMAP-2', cex = .4)
legend("topright", legend = names(grad), col = grad, pch = 16, bty = "n", cex = 2/3)

library('slingshot')
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$ident,
                 start.clus = 'HFSC_1', approx_points = 150)

plot(reducedDims(sce)$UMAP, col = grad[cut(scores$scaled_scores, breaks = 10)],
     asp = 2, pch = 16, xlab = 'UMAP-1', ylab = 'UMAP-2', cex = .4)
lines(SlingshotDataSet(sce), lwd =2, col = 'black')

ks.test(slingPseudotime(sce)[colData(sce)$Assignment == "CMO301", 1],
        slingPseudotime(sce)[colData(sce)$Assignment == "CMO302", 1])

### tradeSeq
"BiocManager::install("tradeSeq")
library(tradeSeq)
set.seed(3)
icMat <- evaluateK(counts = as.matrix(assays(sce)$counts),
                   pseudotime = colData(sce)$slingshot$pseudotime,
                   cellWeights = colData(sce)$slingshot$cellWeights.V1,
                   conditions = factor(colData(sce)$Assignment),
                   nGenes = 300,
                   k = 3:7)

set.seed(3)
sce <- fitGAM(sce, conditions = factor(colData(sce)$Assignment),
              nknots = 5)
mean(rowData(sce)$tradeSeq$converged)

assocRes <- rowData(sce)$assocRes
mockGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionMock, "fdr") <= 0.05)
]
tgfbGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionTGFB, "fdr") <= 0.05)
]

length(mockGenes)"

