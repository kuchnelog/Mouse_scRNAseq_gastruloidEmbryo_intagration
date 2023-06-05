library(SeuratData)
library(ggplot2)
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(kableExtra)
library(stringr)
library(DoubletFinder)
library(cowplot)
source("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/02_scripts/extract_topn.R")
source("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/02_scripts/add_dominantMetadata.R")
source("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/02_scripts/00_generalDeps.R")

# Import Data ------------------------------------------------------------------
em.sing <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/50_rdsObjects/em_sing_process.rds")
gs.sing <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/50_rdsObjects/gs_sing_process.rds")
#gs_em.sing <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/50_rdsObjects/gs_em_sing_cell_idents.rds") 
#gs_em.sing <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/50_rdsObjects/gs_em.rds") 
View(em.sing@meta.data)
View(gs.sing@meta.data)
# Save data path ---------------------------------------------------------------
rdsObject <- "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/50_rdsObjects/"
## table -----------------------------------------------------------------------
table.path <- "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/53_OutPutTable/"
# Genes of interest ------------------------------------------------------------
gene_of_interest <- read.csv("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/genesOfInterest.tsv", header = FALSE,col.names = c("Genes of interest"))
# Parameters -------------------------------------------------------------------
colors.table <- read.table(file="/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/InputTables/ClusterColors.tsv", sep="\t", header=T, comment.char="", as.is=T)
colors.celltype <- setNames(colors.table$blind_friendly[!is.na(colors.table$transferred_identity)], colors.table$transferred_identity[!is.na(colors.table$transferred_identity)])
#  Normalization
norm_method <- "LogNormalize"
# Identification of highly variable features (feature selection)
var_feat_method <- "vst"
# Clustering
Leiden <- 4 
res_var <- 1
general.seed <-42
# Genes of interest
gene_of_interest <- read.csv("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/genesOfInterest.tsv")

# Fusionner les deux objets Seurat--------------------------------------------
features <- SelectIntegrationFeatures(object.list =list(gs.sing, em.sing))
anchorset <- FindIntegrationAnchors(object.list = list(gs.sing, em.sing), dims = 1:20,  anchor.features = features)
gs_em.sing <-  IntegrateData(anchorset = anchorset)

# Identification of highly variable features (feature selection) -------------
gs_em.sing <- FindVariableFeatures(gs_em.sing, selection.method = var_feat_method, nfeatures = 2000)
# Scaling the data -----------------------------------------------------------
all.genes <- rownames(gs_em.sing) 
gs_em.sing <- ScaleData(gs_em.sing, features = all.genes)
# Perform linear dimensional reduction ---------------------------------------
gs_em.sing <- RunPCA(gs_em.sing, features = VariableFeatures(object = gs_em.sing), npcs = 20, nfeatures.print = 10, verbose=TRUE, seed.use=general.seed) # seed.use=NULL
gs_em.sing <- RunUMAP(gs_em.sing, reduction = "pca", dims = 1:20)
gs_em.sing <- FindNeighbors(gs_em.sing, dims = 1:20, verbose=FALSE)
gs_em.sing <- FindClusters(gs_em.sing, resolution = res_var, algorithm = Leiden, random.seed = general.seed)
DimPlot(gs_em.sing, group.by = "orig.ident")
DimPlot(gs_em.sing, reduction = "umap", label = TRUE) + NoLegend()
View(gs_em.sing)

## Supresion des donnes de metadata --------------------------------------------
gs_em.sing$RNA_snn_res.1 <- NULL
gs_em.sing$RNA_snn_res.1.5 <- NULL
gs_em.sing$RNA_snn_res.2 <- NULL
gs_em.sing$RNA_snn_res.3 <- NULL
gs_em.sing$pANN_0.25_0.3_485 <- NULL
gs_em.sing$pANN_0.25_0.1_562 <- NULL
gs_em.sing$DF.classifications_0.25_0.3_485 <- NULL
gs_em.sing$DF.classifications_0.25_0.1_562 <- NULL

# Cell identities --------------------------------------------------------------
atlas.subset <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/50_rdsObjects/atlas_preprocessed2.rds")
anchors <- FindTransferAnchors(reference = atlas.subset, query = gs_em.sing,
                               dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = atlas.subset$celltype,
                            dims = 1:20)
gs_em.sing <- AddMetaData(gs_em.sing, metadata = predictions[,"predicted.id"], col.name = "celltype_integr")
Idents(gs_em.sing) <- gs_em.sing@meta.data$celltype_integr
DimPlot(gs_em.sing)
DimPlot(gs_em.sing,
        pt.size = 1,
        repel = TRUE,
        label = TRUE,
        label.size = 6,
        cols = colors.celltype[levels(Idents(gs_em.sing))]) +
  ggtitle( "Cell identities after Integration") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        line = element_blank()) +
  NoLegend()
saveRDS(gs_em.sing, file = paste0(rdsObject, "gs_em_sing_process.rds"))
rm(list=c("atlas.subset", "anchors", "anchorset", "em.sing", "gs.sing", "predictions"))
gc()





# Heatmap ----------------------------------------------------------------------
DefaultAssay(gs_em.sing) <-"RNA"
#DefaultAssay(gs_em.sing) <-"integrated"
Idents(gs_em.sing)<- gs_em.sing@meta.data$integrated_snn_res.1
gs_em.sing.markers <- FindAllMarkers(gs_em.sing, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gs_em.sing.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
test<-ScaleData(gs_em.sing, verbose = FALSE)
DoHeatmap(test, features = top10$gene) + NoLegend() 

saveRDS(gs_em.sing, file = paste0(rdsObject, "gs_em.rds"))

# 
 # Find markers -----------------------------------------------------------------
DefaultAssay(gs_em.sing) <-"RNA"
 Idents(gs_em.sing) <- gs_em.sing@meta.data$integrated_snn_res.1
gs_em.sing.markers <- FindAllMarkers(gs_em.sing, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
 top10 <- gs_em.sing.markers %>%
   group_by(cluster) %>%
   slice_max(n = 5, order_by = avg_log2FC)
DoHeatmap(gs_em.sing, features = top10$gene) + NoLegend() 


# Idents(gs_em.sing) <- gs_em.sing@meta.data$integrated_snn_res.1
# gs_em.sing.cardio <- subset(gs_em.sing, ident = c(13,14,19,27)) 
# gs_em.sing.markers.cardio <- FindAllMarkers(gs_em.sing.cardio, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# top10 <- gs_em.sing.markers.cardio %>%
#   group_by(cluster) %>%
#   top_n(n = 5, wt = avg_log2FC) 
# DoHeatmap(gs_em.sing.cardio, features = top10$gene) + NoLegend()

Idents(gs_em.sing) <- gs_em.sing@meta.data$integrated_snn_res.1
gs_em.sing.cardio <- subset(test, ident = c(15,5,25)) 
gs_em.sing.markers.cardio <- FindAllMarkers(gs_em.sing.cardio, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- gs_em.sing.markers.cardio %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) 
DoHeatmap(gs_em.sing.cardio, features = top10$gene) + NoLegend()


#markers.sign <- markers[markers$p_val_adj < 0.001 & markers$avg_logFC > log(1.5), ]
gs_em.sing.markers.cardio <- gs_em.sing.markers.cardio[order(-gs_em.sing.markers.cardio$avg_log2FC), ]
gs_em.sing.markers.cardio <- gs_em.sing.markers.cardio[order(gs_em.sing.markers.cardio$cluster), ]
gs_em.sing.markers.cardio <- gs_em.sing.markers.cardio[ , c(6, 7, 2:4, 1, 5)]
View(gs_em.sing.markers.cardio)
gs_em.sing.markers.cardio <- extract_topn(em.sing.cardio.markers)
View(gs_em.sing.markers.cardio)
names_all <- rownames(gs_em.sing.markers.cardio)
names5 <- names_all[1:6]
names15 <- c("Ptn", "Actb", "Actg1", "Vsnl1", "Cald1", "Acta2")
names25 <- c("Mylpf", "Tnnc2", "Tpm2", "Myh3", "Myog", "Tnnt1")

violinFeatureByCluster( names14, em.sing)
VlnPlot(gs_em.sing, names5, ident = c(5,15,25))
VlnPlot(gs_em.sing, names15, ident = c(5,15,25))
VlnPlot(gs_em.sing, names25, ident = c(5,15,25))



# PLOTS ------------------------------------------------------------------------
DimPlot(gs_em.sing, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(gs_em.sing, reduction = "umap", label = TRUE, repel = TRUE,  group.by = "orig.ident") + NoLegend()
DimPlot(object = gs_em.sing, reduction = "umap", group.by = "celltype_DF", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
FeaturePlot(gs_em.sing, features = c("Mylpf","Myh3", "Myh6"),  reduction = "umap", split.by = 'orig.ident', cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))
FeaturePlot(gs_em.sing, features = c("Myh7", "Myl2", "Myl3", "Myl7"),  reduction = "umap", split.by = 'orig.ident', cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))


VlnPlot(gs_em.sing, features="Ftl1-ps1", sort = TRUE, split.by = 'orig.ident') + NoLegend()
VlnPlot(gs_em.sing, features="Ftl1", sort = TRUE, split.by = 'orig.ident') + NoLegend()
VlnPlot(gs_em.sing, features="GFP", sort = TRUE) + NoLegend()
VlnPlot(gs_em.sing, features="Cre", sort = TRUE) + NoLegend()

DefaultAssay(gs_em.sing) <-"RNA"
Idents(gs_em.sing) <- gs_em.sing@meta.data$RNA_snn_res.1
violinFeatureByCluster("GFP", gs_em.sing)
violinFeatureByCluster("Cre", gs_em.sing)

# # Identify conserved cell type markers ---------------------------------------

theme_set(theme_cowplot())
cell25 <- subset(gs_em.sing, idents = 25)
Idents(cell25) <- "orig.ident"
avg.cell25 <- as.data.frame(log1p(AverageExpression(cell25, group.by = "orig.ident",verbose = FALSE)$RNA))
avg.cell25$gene <- rownames(avg.cell25)
genes.to.label = c("Malat1", "Hbb-y", "Ftl1-ps1",  "Ftl1", "Hbb-y", "H19", "Rps2", "Ptma", "Myod1", "Myl7", "Myh8", "Hmgb2", "Xist", "Tagln", "Cdkn1c", "Col2a1", "Gm10260", "Nppb", "Mylpf")
p1 <-  ggplot(avg.cell25, aes(embryE9, gastD11)) + geom_point()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE) + ggtitle("Average Expression: Cluster 25")
p1
View(avg.cell25)

theme_set(theme_cowplot())
cell5 <- subset(gs_em.sing, idents = 5)
Idents(cell5) <- "orig.ident"
avg.cell5 <- as.data.frame(log1p(AverageExpression(cell5, group.by = "orig.ident",verbose = FALSE)$RNA))
avg.cell5$gene <- rownames(avg.cell5)
genes.to.label = c("Malat1", "Hbb-y", "Ftl1-ps1",  "Ftl1", "Hbb-y", "H19", "Rps2", "Ptma", "Hba-x", "Myl7", "Cox6a2", "Hmgb2", "Xist", "Tagln", "Acta1", "Col2a1", "Gm10260", "Nppb", "Acta2")
p1 <-  ggplot(avg.cell5, aes(embryE9, gastD11)) + geom_point()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE) + ggtitle("Average Expression: Cluster 5")
p1
View(avg.cell5)

theme_set(theme_cowplot())
cell15 <- subset(gs_em.sing, idents = 15)
Idents(cell15) <- "orig.ident"
avg.cell15 <- as.data.frame(log1p(AverageExpression(cell15, group.by = "orig.ident",verbose = FALSE)$RNA))
avg.cell15$gene <- rownames(avg.cell15)
genes.to.label = c("Malat1", "Hbb-y", "Ftl1-ps1",  "Ftl1", "Hbb-y", "H19", "Rps2", "Ptma", "Hba-x", "Myl7", "Cox6a2", "Hmgb2", "Xist",  "Gm10260", "Nppb")
p1 <-  ggplot(avg.cell15, aes(embryE9, gastD11)) + geom_point()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE) + ggtitle("Average Expression: Cluster 15")
p1
View(avg.cell15)



View(table1)
table <- table(gs_em.sing@meta.data$celltype_integr, gs_em.sing@meta.data$orig.ident)
table1 <- table(gs_em.sing@meta.data$integrated_snn_res.1, gs_em.sing@meta.data$orig.ident)
View(gs_em.sing@meta.data)
write.table(table1, file = paste0(table.path, "numb_of_cell_clusteurs.csv"))
write.table(table, file = paste0(table.path, "numb_of_cell_origins.csv"))

DefaultAssay(gs_em.sing) <-"integrated"
FeaturePlot(gs_em.sing, features = c("Ftl1","Ftl1-ps1"), split.by = "orig.ident",  reduction = "umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))
View(gs_em.sing@meta.data)          
