library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(kableExtra)
library(stringr)
library(DoubletFinder)
source("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/02_scripts/extract_topn.R")
source("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/02_scripts/add_dominantMetadata.R")
source("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/02_scripts/00_generalDeps.R")

# 1. Load the dataset ----------------------------------------------------------
em.data <- Read10X(data.dir = "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/00_rawData/embryoE9_reAligned",  gene = 2)
em <- CreateSeuratObject(counts = em.data, project = "embryE9")
em.brut <- CreateSeuratObject(counts = em.data, project = "embryE9")

# 2. Parametres ----------------------------------------------------------------
## rds -------------------------------------------------------------------------
rdsObject <- "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/50_rdsObjects/"
## table -----------------------------------------------------------------------
table.path <- "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/53_OutPutTable/"
## QC --------------------------------------------------------------------------
mito <- 5
ribo <- 15
min_counts <- 5000
max_counts <- 90000
min_feats <- 150
min.cells <- 3
# Normalization
norm_method <- "LogNormalize"
# Identification of highly variable features (feature selection)
var_feat_method <- "vst"
## Clustering ------------------------------------------------------------------
Leiden <- 4 
res_var <- c(1, 1.5, 2, 3)
general.seed <- 42
## FindDoublets ----------------------------------------------------------------
dblt.rate <-8
colors.table <- read.table(file="/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/InputTables/ClusterColors.tsv", sep="\t", header=T, comment.char="", as.is=T)
colors.celltype <- setNames(colors.table$blind_friendly[!is.na(colors.table$transferred_identity)], colors.table$transferred_identity[!is.na(colors.table$transferred_identity)])

## Genes of interest -----------------------------------------------------------
gene_of_interest <- read.csv("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/genesOfInterest.tsv", header = TRUE)

# 3. QC ------------------------------------------------------------------------
df <- data.frame(rowSums(em@assays$RNA@counts != 0))
df$features <- rownames(df)
colnames(df) <- c('Nbr_of_cells', 'features')
rownames(df) <- NULL
df_non0 <- df[df$Nbr_of_cells > 0, ]
em[["percent.mt"]] <- PercentageFeatureSet(em, pattern = "^mt-")
em[["percent.ribo"]] <- PercentageFeatureSet(em, pattern="^Rp[sl]")  
em <- subset(em, features=which(df$Nbr_of_cells > min.cells), subset = nCount_RNA >= min_counts & nCount_RNA <= max_counts
             & nFeature_RNA >= min_feats & percent.mt <= mito & percent.ribo > ribo)

##  Visualize QC metrics as a violin plot --------------------------------------
VlnPlot(em, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),  ncol = 4)
## FeatureScatter --------------------------------------------------------------
plot1<- FeatureScatter(em, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +geom_smooth (method= 'lm')
plot2<- FeatureScatter(em, feature1 = "nCount_RNA", feature2 = "percent.mt") +geom_smooth (method= 'lm')
plot3<- FeatureScatter(em, feature1 = "nCount_RNA", feature2 = "percent.ribo") +geom_smooth (method= 'lm')
plot4<- FeatureScatter(em, feature1 = "percent.mt", feature2 = "percent.ribo") +geom_smooth (method= 'lm')
plot1 + plot2 + plot3 + plot4

# 4 Normalization --------------------------------------------------------------
em <- NormalizeData(em, normalization.method = norm_method, scale.factor = 10000)

# 5 Identification of highly variable features (feature selection) -------------
em <- FindVariableFeatures(em, selection.method = var_feat_method, nfeatures = 2000)
## Identify the 10 most highly variable genes ----------------------------------
top10 <- head(VariableFeatures(em), 10)
## plot variable features with and without labels-------------------------------
plot1 <- VariableFeaturePlot(em)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# 6. Scaling the data -----------------------------------------------------------
all.genes <- rownames(em) 
em <- ScaleData(em, features = all.genes)
# 7. Perform linear dimensional reduction ---------------------------------------
em <- RunPCA(em, features = VariableFeatures(object = em), npcs = 20, nfeatures.print = 10, verbose=TRUE, seed.use=general.seed)
ElbowPlot(em, ndims = 20)

em <- FindNeighbors(em, dims = 1:20, verbose=FALSE)
em <- FindClusters(em, resolution = res_var, algorithm = Leiden, random.seed = general.seed) # algorithm = algo.cluster, random.seed = general.seed, verbose=FALSE
em <- RunUMAP(em, dims = 1:20, seed.use = general.seed)
DimPlot(em, group.by = "RNA_snn_res.1" , reduction = "umap", label = TRUE) + NoLegend()
#8.  Finding differentially expressed features (cluster biomarkers)--------------
em.markers <- FindAllMarkers(em, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, group.by = "RNA_snn_res.1")
em.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(em, features = top10$gene, group.by = "RNA_snn_res.1") + NoLegend()

# 9. Doubets Finder ------------------------------------------------------------
atlas.subset <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/50_rdsObjects/atlas_preprocessed2.rds")
# atlas.subset <- readRDS(paste0(rdsObject,"atlas_preprocessed2.rds")) # reuse the variable you created at the beginning
anchors <- FindTransferAnchors(reference = atlas.subset, query = em,
                               dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = atlas.subset$celltype,
                            dims = 1:20)
em <- AddMetaData(em, metadata = predictions[,"predicted.id"], col.name = "celltype_DF")
Idents(em) <- em@meta.data$celltype_DF

DimPlot(em,
        pt.size = 1,
        repel = TRUE,
        label = TRUE,
        label.size = 6,
        cols = colors.celltype[levels(Idents(em))]) +
  ggtitle( " Cell identities to perform DoubletFinder") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        line = element_blank()) +
  NoLegend()
#
View(em@meta.data)
## pK identifikation -----------------------------------------------------------
sweep.res <- paramSweep_v3(em, PCs = 1:20) # as estimated from PC elbowPlot
sweep.stats_em <- summarizeSweep(sweep.res, GT = FALSE)
Sys.sleep(0.5)
bcmvn_em <- find.pK(sweep.stats_em)

ggplot(bcmvn_em, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_em %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 

pK <- as.numeric(as.character(pK[[1]]))
## Homotypic Doublets Proportion Estimate --------------------------------------
annotations <- em@meta.data$celltype_DF
homotypic.prop <- modelHomotypic(annotations)
nDoublets <- round(ncol(em)*dblt.rate/100) # poi
nDoublets_nonhomo <- round(nDoublets*(1-homotypic.prop)) #poi.adj

## Run doubletFinder: creating artificial doublets------------------------------
em <- doubletFinder_v3(em,
                       PCs = 1:20,
                       pN = 0.25,
                       pK = pK,
                       nExp = nDoublets_nonhomo)

## Plot the singlets, doublets and the count of UMIs in each cells--------------
col_dblts <- grep("DF.classifications", colnames(em@meta.data), value=TRUE)
Idents(em) <- col_dblts

cellsData <- data.frame(em@reductions[["umap"]]@cell.embeddings, em@meta.data[,col_dblts])
colnames(cellsData) <- c(colnames(em@reductions[["umap"]]@cell.embeddings), "col_dblts")
p1 <- ggplot(cellsData[c('UMAP_1', 'UMAP_2')], # Omit the column used for facetting to get all points repeated in all facets
             aes( x = UMAP_1,
                  y = UMAP_2)) +
  geom_point( alpha = .4,
              size  = .75,
              color = "grey") +
  geom_point( data = cellsData, # Now provide data including column for facetting
              aes(color = col_dblts),
              alpha = .5,
              size  = 1.2) +
  facet_wrap(facets = vars(factor(col_dblts, levels = c("Singlet", "Doublet"))),
             ncol = 2) +
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        line = element_blank())
p1
# 
cellsData <- data.frame(em@reductions[["umap"]]@cell.embeddings, em@meta.data$nCount_RNA)
colnames(cellsData) <- c(colnames(em@reductions[["umap"]]@cell.embeddings), "nCount_RNA")
Db <- table(em@meta.data[col_dblts])["Doublet"]
Sing <- table(em@meta.data[col_dblts])["Singlet"]
saveRDS(em, file = paste0(rdsObject, "em.rds"))
#em <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/50_rdsObjects/em.rds")
## Subset to singlets ----------------------------------------------------------
em.sing <- subset(em, idents='Singlet')
saveRDS(em.sing, file = paste0(rdsObject, "em_singlets.rds"))

# 10. Preprocessing workflow (after doublets removal) -------------------------
em.sing <- FindVariableFeatures(em.sing, selection.method = var_feat_method, nfeatures=2000, verbose=FALSE)
#var.feats.dataset <- VariableFeatures(em.sing)
em.sing <- ScaleData(em.sing, features=rownames(em.sing), do.scale=FALSE, verbose=FALSE)
em.sing <- RunPCA(em.sing, npcs = 20, nfeatures.print = 10, seed.use = general.seed, verbose=TRUE)
# 11. Clustering ----------------------------------------------------------------
em.sing <- FindNeighbors(em.sing, dims = 1:20, verbose=FALSE)
em.sing <- FindClusters(em.sing, resolution = res_var, algorithm = Leiden, random.seed = general.seed) # algorithm = algo.cluster, random.seed = general.seed, verbose=FALSE
## Run non-linear dimensional reduction (UMAP) ---------------------------------
em.sing <- RunUMAP(em.sing, dims = 1:20, seed.use = general.seed)
DimPlot(em.sing, group.by = "RNA_snn_res.1" , reduction = "umap", label = TRUE) + NoLegend()

# 12. Finding differentially expressed features (cluster biomarkers)-------------
em.sing.markers <- FindAllMarkers(em.sing, group.by = "RNA_snn_res.1", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- em.sing.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) 
DoHeatmap(em.sing, features = top10$gene, group.by = "RNA_snn_res.1") + NoLegend()

em.sing.markers <- FindAllMarkers(em.sing, group.by = "RNA_snn_res.3", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- em.sing.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) 
DoHeatmap(em.sing, features = top10$gene, group.by = "RNA_snn_res.3") + NoLegend()
DimPlot(em.sing, group.by = "RNA_snn_res.3" , reduction = "umap", label = TRUE) + NoLegend()

Idents(em.sing) <- em.sing@meta.data$RNA_snn_res.3
em.sing.cardio <- subset(em.sing, ident = c(15,44,41,45)) 
em.sing.cardio.markers <- FindAllMarkers(em.sing.cardio, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- em.sing.cardio.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) 
DoHeatmap(em.sing.cardio, features = top10$gene) + NoLegend()

#markers.sign <- markers[markers$p_val_adj < 0.001 & markers$avg_logFC > log(1.5), ]
em.sing.cardio.markers.sign <- em.sing.cardio.markers[order(-em.sing.cardio.markers$avg_log2FC), ]
em.sing.cardio.markers.sign <- em.sing.cardio.markers.sign[order(em.sing.cardio.markers.sign$cluster), ]
em.sing.cardio.markers.sign <- em.sing.cardio.markers.sign[ , c(6, 7, 2:4, 1, 5)]
View(em.sing.cardio.markers.sign)

em.sing.cardio.markers <- extract_topn(em.sing.cardio.markers)
names15 <- rownames(em.sing.cardio.markers.sign)
names15 <- names15[1:6]
names41 <- c("Tnnc2", "Meg3", "Itm2a", "Myh3", "Cdkn1c", "Tpm2")
names44 <- c("Myh6", "Crip1", "Tnni3", "Stard10", "5430431A17Rik", "Atp2a2")
names45 <- c("Myl2", "Cdk1", "Cited1", "Ube2c", "Ccnb1", "Cks2")
violinFeatureByCluster( names14, em.sing)
VlnPlot(em.sing, names15, ident = c(15,41,44,45))
VlnPlot(em.sing, names41, ident = c(15,41,44,45))
VlnPlot(em.sing, names44, ident = c(15,41,44,45))
VlnPlot(em.sing, names45, ident = c(15,41,44,45))

## Table de marquers -----------------------------------------------------------
write.table(em.sing.cardio.markers.sign, file = paste0(table.path, "marquers_cl_14_41_44_45_em_sign.csv"))
write.table(cluster14_35.markers, file = paste0(table.path, "marquers_cl_14_35_em.csv"))
write.table(cluster14_45.markers, file = paste0(table.path, "marquers_cl_14_45_em.csv"))
# VlnPlot(em.sing, cluster14_35.markers.14names[1:3], idents = c(14, 35))
# VlnPlot(em.sing, cluster14_35.markers.35names[1:9], idents = c(14, 35))
# VlnPlot(em.sing, cluster14_45.markers.14names[1:9],  idents = c(14, 45))
# VlnPlot(em.sing, cluster14_45.markers.45names[1:2],  idents = c(14, 45))

# 13. Cell identities after perform DoubletFinder--------------------------------
anchors <- FindTransferAnchors(reference = atlas.subset, query = em.sing,
                               dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = atlas.subset$celltype,
                            dims = 1:20)
em.sing <- AddMetaData(em.sing, metadata = predictions[,"predicted.id"], col.name = "celltype_sing")
Idents(em.sing) <- em.sing@meta.data$celltype_sing

DimPlot(em.sing,
        pt.size = 1,
        repel = TRUE,
        label = TRUE,
        label.size = 6,
        cols = colors.celltype[levels(Idents(em.sing))]) +
  ggtitle( "Cell identities") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        line = element_blank()) +
  NoLegend()

View(atlas.subset)
saveRDS(em.sing, file = paste0(rdsObject,"em_sing_process.rds"))
## Celltype identity
Idents(em.sing) <- factor(em.sing@meta.data$celltype, levels = sort(levels(as.factor(em.sing@meta.data$celltype))))

## Plots -----------------------------------------------------------------------

FeaturePlot(em.sing, features = c("Ftl1","Ftl1-ps1"),  reduction = "umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))
FeaturePlot(em.sing, features = c("GFP","Cre"),  reduction = "umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))
violinFeatureByCluster("GFP", em.sing)
violinFeatureByCluster("Cre", em.sing)


saveRDS(em.sing, file = paste0(rdsObject,"em_sing_total.rds"))
em.sing <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/50_rdsObjects/em_sing_total.rds")
