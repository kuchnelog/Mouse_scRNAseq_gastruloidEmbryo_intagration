library(Seurat)
library(Matrix)
library(harmony)
library(cowplot)

# Path management
basePath <- "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/seuratAnalysis"
inputs.folder <- file.path(basePath, "inputData")
atlas.folder <- file.path(inputs.folder, "atlas")

# Load data and create Seurat object
atlas.metadata <- read.table(file.path(atlas.folder, "meta.tab") , sep="\t" , header=TRUE )
atlas.genes <- read.csv(file.path(atlas.folder, "genes.tsv"), sep="\t", header = FALSE, as.is = TRUE)
atlas.data <- readMM(file.path(atlas.folder, "raw_counts.mtx"))
atlas.data@Dimnames[[2]] <- atlas.metadata$cell
atlas.data@Dimnames[[1]] <- atlas.genes$V2

rm(atlas.genes)
gc()


atlas <- CreateSeuratObject(counts = atlas.data, project = "Atlas", min.cells=3, min.features = 200)
rm(atlas.data)
gc()

# Reduce metadata table to the list of cell in the atlas
nonFiltered <- intersect(colnames(atlas), atlas.metadata$cell)
atlas.metadata <- atlas.metadata[which(atlas.metadata$cell %in% nonFiltered),]
gc()

# Normalize and add other metadata
atlas <- NormalizeData(atlas, normalization.method = "LogNormalize", verbose=FALSE)
atlas <- AddMetaData(object = atlas, metadata = atlas.metadata$stage, col.name = "day")
atlas <- AddMetaData(object = atlas, metadata = atlas.metadata$celltype, col.name = "celltype")
atlas <- AddMetaData(object = atlas, metadata = atlas.metadata$sequencing.batch, col.name = "replicate")
atlas <- AddMetaData(object = atlas, metadata = paste(atlas$day, "_batch_0", atlas$replicate, sep=""), col.name = "dataset")
atlas$model <- "Embryos"
rm(atlas.metadata, nonFiltered)
gc()

# preprocess the atlas
atlas.subset <- FindVariableFeatures(atlas, nfeatures=2000, selection.method = "vst", verbose=FALSE)
atlas.subset <- ScaleData(atlas.subset, do.scale = FALSE, verbose = FALSE)
atlas.subset <- RunPCA(atlas.subset, npcs = 50, verbose = FALSE)
gc()

p1 <- DimPlot(object = atlas.subset, reduction = "pca", pt.size = .1, group.by = "day")
p2 <- VlnPlot(object = atlas.subset, features = "PC_1", group.by = "day", pt.size = .1)
plot_grid(p1,p2)

# Save the preprocessed atlas
#saveRDS(atlas.subset, file.path(atlas.folder, paste0("atlas_preprocessed2.rds")))
rdsObject <- "/mnt/DATA_4TB/projects/gastruloids_sc_Lescroart/analysis/embryE9gastD10/50_rdsObjects/"
saveRDS(atlas.subset, file = paste0(rdsObject, "atlas_preprocessed2.rds"))

# Create Umap of the atlas
top.pcs <- 30
general.seed <- 17
atlas.subset <- RunUMAP(atlas.subset, dims = 1:top.pcs, seed.use = general.seed, verbose = FALSE)

saveRDS(atlas.subset, file.path(atlas.folder, paste0("atlas_preprocessed_UMAP.rds")))