#Load Libraries 
library(Seurat)
library(readr)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)

#####---- GSE129730 Seurat Object----####

# Path and file name of downloaded data - Change this to whatever your path is
count_file_MB_GSE129730 <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE129730_mouse_primary_MB_Drop_2022_05_27\\GSE129730_scRNAseq_2022_05_27\\GSE129730_RAW\\GSM3720938_WT03.dge.txt.gz"

# Read in data as a dataframe
counts_MB_GSE129730 <- as.data.frame(read_delim(count_file_MB_GSE129730,    #file to read in
                                   delim = '\t')) #Tab delimiter - how columns are separated
## Replace NA's with 0
counts_MB_GSE129730[is.na(counts_MB_GSE129730)] <- 0

# Make gene names Row name
rownames(counts_MB_GSE129730) <- counts_MB_GSE129730$GENE
counts_MB_GSE129730 <- counts_MB_GSE129730[-1]

#Create a seuratobject
seurat_object_GSE129730 <-  CreateSeuratObject(counts_MB_GSE129730, project = "GSE129730")



#####---- GSE175416 Seurat Object ----#####
# Assign file paths to mtx, cells, features to read in files separately

mtx <- "J:\\\\Biostat_Interns\\\\Nobles_Gabrielle\\\\Moffitt_Internship_Data_2022\\\\Internship_scRNAseq_Data_2022\\\\Mouse_scRNAseq_Data\\\\GSE175416_mouse_primary_MB_10x_2022_05_27\\\\GSE175416_scRNA_seq_2022_05_27\\\\GSE175416_RAW\\\\GSM5333101_Primary_Tumor-matrix.mtx.gz"
cells <- "J:\\\\Biostat_Interns\\\\Nobles_Gabrielle\\\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE175416_mouse_primary_MB_10x_2022_05_27\\GSE175416_scRNA_seq_2022_05_27\\GSE175416_RAW\\GSM5333101_Primary_Tumor-barcodes.tsv.gz"
features <- "J:\\\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE175416_mouse_primary_MB_10x_2022_05_27\\GSE175416_scRNA_seq_2022_05_27\\GSE175416_RAW\\GSM5333101_Primary_Tumor-genes.tsv.gz"

# ReadMtx loads in sparse data 
counts_MB_GSE175416 <- ReadMtx(mtx, cells, features)

#Create a seuratobject 
seurat_object_GSE175416 = CreateSeuratObject(counts_MB_GSE175416, 
                                             project = "GSE175416")




######---- MB_Combined QC----####
# Merge the seurat objects 
MB_combined <- merge(seurat_object_GSE129730, y = seurat_object_GSE175416, add.cell.ids = c("MB_GSE129730", "MB_GSE175416"), project = "Combined_MB")

# Normalize the MB_combined dataset 
# MB_combined <- NormalizeData(MB_combined, normalization.method = "CLR")

# Find high variable features of MB_combined
MB_combined <- FindVariableFeatures(MB_combined, selection.method = "vst", 
                                      nfeatures = 2000)

# Identify Top 10 HVG of MB_combined 
top10 <- head(VariableFeatures(MB_combined), 10)
top10

# ScaleData converts normalized gene expression to Z-score of MB_combined 
all.genes <- rownames(MB_combined)
MB_combined <- ScaleData(MB_combined, features = all.genes)


##----MB_combined Dimensional Reduction----##
# Run PCA for linear dimensional reduction 
MB_combined <- RunPCA(MB_combined,
                        features = VariableFeatures(object = MB_combined), 
                      seed.use = 42)


####----MB_combined Cluster the Cells----####
# Compute nearest neighbor grqph and SNN: shared nearest neighbor 
MB_combined <- FindNeighbors(MB_combined, dims = 1:10)

# Find clusters 
MB_combined <- FindClusters(MB_combined, resolution = 0.5)


####----MB_combined tSNE & UMAP----####
# Generate a UMAP for visualization 
MB_combined <- RunTSNE (MB_combined, dims = 1:10, check_duplicates = FALSE, 
                        seed.use = 1)

MB_combined <- RunUMAP(MB_combined, dims = 1:10, check_duplicates = FALSE,
                       seed.use = 42) 


# DimPlot use UMAP by default, Seurat clusters as identity  
DimPlot(MB_combined, reduction = "tsne",group.by = "orig.ident", seed = 1)

DimPlot(MB_combined, reduction = "umap",group.by = "orig.ident", seed = 1)

?DimPlot()





####----MB.Normalized merge seurat object ----####
# Normzlize GSE129730 and GSE175416 seurat objects 
seurat_object_GSE129730_normalized <- NormalizeData(seurat_object_GSE129730,
                                                    normalization.method = "CLR")

seurat_object_GSE175416_normalized <- NormalizeData(seurat_object_GSE175416,
                                                    normalization.method = "CLR")

# Merged normalized data 
MB.normalized <- merge(seurat_object_GSE129730_normalized,
                       y = seurat_object_GSE175416_normalized,
                         add.cell.ids = c("GSE129730_normalized", "GSE175416_normalized"),
                       project = "Normalized ",merge.data = TRUE)


# Find high variable features 

MB.normalized <- FindVariableFeatures(MB.normalized, selection.method = "vst", 
                                      nfeatures = 2000)

# Identify Top 10 HVG
top10 <- head(VariableFeatures(MB.normalized), 10)
top10

# ScaleData converts normalized gene expression to Z-score 
all.genes <- rownames(MB.normalized)
MB.normalized <- ScaleData(MB.normalized, features = all.genes)


####---- MB.Normalized Dimensional Reduction----####
# Run PCA for linear dimensional reduction 
MB.normalized <- RunPCA(MB.normalized,
                        features = VariableFeatures(object = MB.normalized), 
                        seed.use = 42)

####---- MB.Normalized Cluster the Cells----####
# Compute nearest neighbor grqph and SNN: shared nearest neighbor 
MB.normalized <- FindNeighbors(MB.normalized, dims = 1:10)

# Find clusters 
MB.normalized <- FindClusters(MB.normalized, resolution = 0.5)


####----tSNE----####
# Generate a TSNE and Umap for visualization 
MB.normalized <- RunTSNE (MB.normalized, dims = 1:10, check_duplicates = FALSE,
                          seed.use = 1) 

MB.normalized <- RunUMAP(MB.normalized, dims = 1:10, check_duplicates = FALSE,
                         seed.use = 42) 


# DimPlot use tsne by default, Seurat clusters as identity  
DimPlot(MB.normalized, reduction = "tsne", group.by = "orig.ident",seed = 1) 

DimPlot(MB.normalized,label= T, repel = F, label.size = 4, group.by = "seurat_clusters", seed = 1)



#####---- Cell Annotation ----####

#Select the functions 
# hpca.ref <- celldex:: HumanPrimaryCellAtlasData()
monaco.ref <- celldex::MonacoImmuneData()
mouseRNA.ref <- celldex::MouseRNAseqData()


# Convert Seurat object to a single cell experiment for convince 
sce.normalized <- as.SingleCellExperiment(DietSeurat(MB.normalized))


sce.not.normalized <- as.SingleCellExperiment((DietSeurat(MB_combined)))

#mouseRNAseq data sce.normalized 
mouse.main <- SingleR(test = sce.normalized,assay.type.test = 1,ref = mouseRNA.ref,
                       labels = mouseRNA.ref$label.main)
mouse.fine <- SingleR(test = sce.normalized,assay.type.test = 1,ref = mouseRNA.ref,
                       labels = mouseRNA.ref$label.fine)

#mouseRNAseq data sce.not.normalzied
mouse.main <- SingleR(test = sce.not.normalized,assay.type.test = 1,ref = mouseRNA.ref,
                      labels = mouseRNA.ref$label.main)
mouse.fine <- SingleR(test = sce.not.normalized,assay.type.test = 1,ref = mouseRNA.ref,
                      labels = mouseRNA.ref$label.fine)


# mouseRNA MB.normalized 
MB.normalized@meta.data$mouse.main <- mouse.main$pruned.labels
MB.normalized@meta.data$mouse.fine <- mouse.fine$pruned.labels

# mouseRNA MB_combined 
MB_combined@meta.data$mouse.main <- mouse.main$pruned.labels
MB_combined@meta.data$mouse.fine <- mouse.fine$pruned.labels


# mouseRNA MB.normalized
MB.normalized<- SetIdent(MB.normalized, value = "mouse.fine")

# mouseRNA MB_combind
MB_combined <- SetIdent(MB_combined, value = "mouse.fine")

## Tsne MB.normalized
DimPlot(MB.normalized, label = T , repel = F,
        label.size = 2, reduction = "tsne" ,seed = 1) 
  

DimPlot(MB.normalized, label = T, repel = F, 
        label.size = 3, reduction = "umap", seed = 1)

## tsne MB_combined
DimPlot(MB_combined, label = T , repel = F,
        label.size = 3, reduction = "tsne", seed = 1) 

DimPlot(MB_combined, label = T, repel = F, 
        label.size = 3, reduction = "umap", seed = 1) 



####----Differential Gene expression----####

all.markers <- FindAllMarkers(MB.normalized, only.pos = T, min.pct = 0.5, 
                              logfc.threshold = 0.5)
# To find the 
top5.markers<- as.data.frame(all.markers %>% group_by(cluster) %>% 
                                top_n(n = 5, wt = avg_log2FC))
top5.markers


Microglia.plot <- FeaturePlot(MB.normalized, 
                              reduction = "umap", 
                              features = c("Itgam","P2ry12", "Tmem119", "Cx3cr1"),
                              label = FALSE, 
                              order = TRUE,
                              
) 


Microglia.plot
FeaturePlot()


Astrocyte.plot <- FeaturePlot(MB.normalized, 
                              reduction = "umap", 
                              features = c("S100b","Slc1a3", "Aldh1l1", "Gfap"),
                              label = FALSE, 
                              order = TRUE,
                              label.size = 1.60
) 

Astrocyte.plot

Fibroblast.plot <- FeaturePlot(MB.normalized, 
                            reduction = "umap", 
                            features = c("Fn1","Tnc"),
                            label = FALSE, 
                            order = TRUE,
                            label.size = 1.60
) 

Fibroblast.plot

Tanycytes.plot <- FeaturePlot(MB.normalized, 
                               reduction = "umap", 
                               features = c("Gpr50","Col23a1"),
                               label = FALSE, 
                               order = TRUE,
                               label.size = 1.60
 ) 

Tanycytes.plot

Oligo.precursor.plot <- FeaturePlot(MB.normalized, 
                              reduction = "umap", 
                              features = c("Lhfpl3","Megf11", "Pcdh15", "Cspg4"),
                              label = FALSE, 
                              order = TRUE,
                              label.size = 1.60
) 

Oligo.precursor.plot

oligodendrocyte.plot <-  FeaturePlot(MB.normalized, 
                                reduction = "umap", 
                                features = c("Mog","Mbp","Mag","Cldn11"),
                                label = FALSE, 
                                order = TRUE,
                                label.size = 1.60
) 

oligodendrocyte.plot

oligodendrocyte.plot2 <- FeaturePlot(MB.normalized,
            reduction = "umap",
            features = c("Olig1","Olig2","Olig3"),
            label = F, 
            order = T, 
            label.size = 1.6)


oligodendrocyte.plot2

immunesystem.plot <-  FeaturePlot(MB.normalized, 
                                  reduction = "umap", 
                                  features = c("Il2ra","Ccr6","Cxcr3","Cd4"),
                                  label = FALSE, 
                                  order = TRUE,
                                  label.size = 1.60
) 

# Neural Progintor Cells 
npc.plot <-  FeaturePlot(MB.normalized, 
                                  reduction = "umap", 
                                  features = c("Cspg4","Rnf5","Sox2","Sox1"),
                                  label = FALSE, 
                                  order = TRUE,
                                  label.size = 1.60
) 

npc.plot

# Neural Stem Cells 
nsc.plot <-  FeaturePlot(MB.normalized, 
                         reduction = "umap", 
                         features = c("Sox9","Prom1","Nes"),
                         label = FALSE, 
                         order = TRUE,
                         label.size = 1.60
) 
nsc.plot



macrophage.plot <-  FeaturePlot(MB.normalized, 
                         reduction = "umap", 
                         features = c("Cd68","Cd163","Cd14"),
                         label = FALSE, 
                         order = TRUE,
                         label.size = 1.60
) 
macrophage.plot




# To write out table for cluster metadata 
write.table(MB.normalized@meta.data, file="C:\\Users\\80023619\\Documents\\MB.normalized.meta.data.txt", sep="\t")



DimPlot(MB.normalized, reduction = "umap", group.by = "seurat_clusters") 


####### remove clusters
# Can also provide multiple idents using the typical R syntax "c()"
sub_obj <- subset(MB.normalized, idents = c("NPCs","aNSCs"), invert = TRUE)


DimPlot(sub_obj, label = T , repel = F,
        label.size = 2, reduction = "tsne" ,group.by = "mouse.fine") 

sub_obj2 <- subset(MB_combined, idents = c("NPCs","aNSCs"), invert = TRUE)

DimPlot(sub_obj2, label = T , repel = F,
        label.size = 1.5, reduction = "tsne" ,group.by = "mouse.fine") 
