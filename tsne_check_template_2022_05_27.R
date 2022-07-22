# Load Library 
library(Seurat)
library(SeuratObject)
library(readr)
library(SeuratDisk)



####----Reading in File----####

# Path and file name of downloaded data - Change this to whatever your path is
count_dir <- "/Users/gabriellenobles/Desktop/scRNA-seq project/GSE70630_OG_processed_data_v2.txt"

#.txt or .csv files 
raw_data <- read.table(file = , header = T, row.names=1, sep="", as.is=T)

# 10X CellRanger .HDF5 format 
h5_data <- Read10X_h5(filename = count_dir,
                       use.names = TRUE,
                       unique.features = TRUE)
# Transform count into a dataframe 
counts <- as.data.frame(read_delim(raw_data,    #file to read in
                                   delim = '\t')) #Tab delimiter - how columns are separated


## If you have sparse 10x files with own data
## Assign file paths to mtx, cells, features to read in files separately

mtx <- "" # matrix  
cells <- "" # barcode
features <- "" # gene
# ReadMtx loads in sparse data 
counts <- ReadMtx(mtx, cells, features)

####---- Replace NA's with 0's ----####
# Replace NA's with 0
counts[is.na(counts)] <- 0
# Make gene names Row name
rownames(counts) <- counts$GENE
counts <- counts[-1]


####---Create a seuratobject----####
seurat_object <- CreateSeuratObject(counts)

# 
# #Calculate Qc metrics for mitochondria
# seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
# 
# head(seurat_object@meta.data, 5) ##Show QC matrix for first 5 cells
# 
# #Visualize QC matrix as a violin plot 
# VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# #Visualize feature-feature relationship between RNA count and mitochondria
# #Visualize feature-feature relationship between RNA count and RNA features(genes)
# plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# 

####----Normalization and Highly variable features----####

# Normalize the data to account for sequencing depth
seurat_object <- NormalizeData(seurat_object, normalization.method = "CLR")

# Find high variable features 
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", 
                                      nfeatures = 2000)

# Identify Top 10 HVG
top10 <- head(VariableFeatures(seurat_object), 10)
top10

# ScaleData converts normalized gene expression to Z-score 
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)


####----Dimensional Reduction----####
# Run PCA for linear dimensional reduction 
seurat_object <- RunPCA(seurat_object,
                        features = VariableFeatures(object = seurat_object),
                        seed.use = 42)

# Examine PCA
print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)

# # Visualizing principal component Loadings with ggplot
# VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca") 
# 
# # Visulazing principal component loadings with heat maps 
# DimHeatmap(seurat_object, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
# 

# DimPlot to visualize all reduced representations 
## PCA, tSNE, UMAP etc. 
# 
# DimPlot(seurat_object, reduction = "pca") 
# 
# ElbowPlot(seurat_object)


####----Cluster the Cells----####
# Compute nearest neighbor grqph and SNN: shared nearest neighbor 
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)

# Find clusters 
seurat_object <- FindClusters(seurat_object, resolution = 0.5)


####----tSNE----####
# Generate a UMAP for visualization 
seurat_object <- RunTSNE (seurat_object, dims = 1:10)

# DimPlot use UMAP by default, Seurat clusters as identity  
DimPlot(seurat_object, reduction = "tsne") 
  
 
  
  