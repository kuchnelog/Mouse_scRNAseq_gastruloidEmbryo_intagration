library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(kableExtra)
library(stringr)
library(DoubletFinder)


# Charger les fichiers----------------------------------------------------------
#cells <- read.csv("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/00_rawData/atlas2/metadata_cells.csv",  header = TRUE)
genes <- read.csv("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/00_rawData/atlas2/metadata_genes.csv", header = TRUE)
#genes_subset <- genes[, c("id", "mgi_symbol")]
named_vector <- setNames(genes$mgi_symbol, genes$id)

#SCE
atlas2 <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/00_rawData/atlas2/embryo_sce.rds")
#atlas2 <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/00_rawData/atlas2/atlas2.rds")
atlas2 <- as.Seurat(atlas2, data = NULL) #Seurat obj

# Parametres--------------------------------------------------------------------
# Normalization
norm_method <- "LogNormalize"
# Identification of highly variable features (feature selection)
var_feat_method <- "mvp"
rdsObject <- "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/00_rawData/atlas2/"

# Peparation de l'atlas avec les donnes reduites et modifies--------------------
Idents(atlas2) <- atlas2$stage
atlas2 <- subset(atlas2, idents = c("E9.0", "E9.5"))
atlas2$leiden <- NULL
atlas2$louvain <- NULL
atlas2$somite_count <-NULL
atlas2$S_score <- NULL
atlas2$G2M_score <- NULL
atlas2$phase <- NULL
colnames(atlas2@meta.data)[3] <- "nFeature_RNA"
colnames(atlas2@meta.data)[2] <- "nCount_RNA"
View(atlas2@meta.data)
#ensemble_ids <- GetAssayData(atlas2, slot = "counts" )@Dimnames[[1]]
#mgi_symbols <- named_vector[ensemble_ids]

str(atlas2)
# Extraire les noms de lignes de la matrice dans le slot "counts" de votre objet Seurat
#rownames_counts <- rownames(atlas2@assays$originalexp@counts)
# Rename the names of the counts matrix with the MGI symbols
#names(atlas2@assays$originalexp@counts@Dimnames[[1]]) <- mgi_symbols
#SetAssayData(atlas2, slot = "counts", assay = "originalexp", data = as.matrix(mgi_symbols))
# Replace the row names of the counts matrix with MGI symbols
#rownames(atlas2@assays$originalexp@counts) <- mgi_symbols
#rownames(atlas2@assays$originalexp@data) <- mgi_symbols
#str(atlas2)
# Update the metadata with the new MGI symbols
#atlas2@meta.data$mgi_symbols <- mgi_symbols

# Normalization ----------------------------------------------------------------
atlas2 <- NormalizeData(atlas2, normalization.method = norm_method, scale.factor = 10000)
# Identification of highly variable features (feature selection) -------------
atlas2 <- FindVariableFeatures(atlas2, selection.method = var_feat_method, nfeatures = 2000)
# Scaling the data -------------------------------------------------------------
all.genes <- rownames(atlas2) 
atlas2 <- ScaleData(atlas2, features = all.genes)
# PCA --------------------------------------------------------------------------
atlas2 <- RunPCA(atlas2, features = VariableFeatures(object = atlas2), npcs = 20, nfeatures.print = 10, verbose=TRUE, seed.use=NULL)
#names(atlas2@reductions$pca@feature.loadings) <- mgi_symbols
#names(atlas2@reductions$pca@feature.loadings[[1]]) <- mgi_symbols
ElbowPlot(atlas2, ndims = 20)
str(atlas2)

saveRDS(atlas2, file = paste0(rdsObject, "atlas2brut.rds"))
#atlas2 <- readRDS("/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/00_rawData/atlas2/atlas2brut.rds")

# Modification des ids des genes -----------------------------------------------
atlas2@assays$originalexp@counts@Dimnames[[1]] <- named_vector
atlas2@assays$originalexp@data@Dimnames[[1]] <- named_vector
atlas2@assays$originalexp@meta.features[[1]] <- named_vector
atlas2@assays$originalexp@var.features <- named_vector ## important, not in rds
atlas2@assays$originalexp@scale.data <- named_matrix
atlas2@reductions$pca@feature.loadings <- named_vector ## il faut le faire - done
# Convert named_vector to matrix
named_matrix <- as.matrix(named_vector)

# Assign named_matrix to feature.loadings in DimReduc object
atlas2@reductions$pca@feature.loadings <- named_matrix



atlas2@reductions$PCA@params$features <- named_vector
atlas2@commands$ScaleData.originalexp@params$features <- named_vector
atlas2@commands$RunPCA.originalexp@params$features <- named_vector
str(atlas2)
# Examine and visualize PCA results a few different ways

print(atlas2[["pca"]], dims = 1:5, nfeatures = 5)
atlas2 <- RunUMAP(atlas2, dims = 1:20, seed.use = NULL)
DimPlot(atlas2, reduction = "umap", label = TRUE) + NoLegend()

# Save object ------------------------------------------------------------------
saveRDS(atlas2, file = paste0(rdsObject, "atlas2.rds"))

# Plots ------------------------------------------------------------------------
id_cell <- rownames(atlas2@meta.data$cell)
Idents(atlas2) <- atlas2@meta.data$celltype_extended_atlas
DimPlot(atlas2,
        pt.size = 1,
        repel = TRUE,
        label = TRUE,
        label.size = 6,
        cols = (Idents(atlas2)) +
  ggtitle( "Cell identities") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        line = element_blank())) + NoLegend()

