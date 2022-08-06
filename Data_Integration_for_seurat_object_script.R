#Load Libraries 
library(Seurat)
library(readr)
library(SingleR)
library(dplyr)
library(celldex)
library(SingleCellExperiment)

#####---- GSE129730 Seurat Object----####

# Path and file name of downloaded data - Change this to whatever your path is
count_file_MB_GSE129730 <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE129730_mouse_primary_MB_Drop_2022_05_27\\GSE129730_scRNAseq_2022_05_27\\GSE129730_RAW\\GSM3720938_WT03.dge.txt.gz"

# Read in data as a dataframe
counts_MB_GSE129730 <- as.data.frame(read_delim(count_file_MB_GSE129730,    #file to read in
                                   delim = '\t')) #Tab delimiter - how columns are separated


#Create a seuratobject
seurat_object_GSE129730 <-  CreateSeuratObject(counts_MB_GSE129730, project = "GSE129730")



#####---- GSE175416 Seurat Object ----#####
# Assign file paths to mtx, cells, features to read in files separately

mtx_1 <- "J:\\\\Biostat_Interns\\\\Nobles_Gabrielle\\\\Moffitt_Internship_Data_2022\\\\Internship_scRNAseq_Data_2022\\\\Mouse_scRNAseq_Data\\\\GSE175416_mouse_primary_MB_10x_2022_05_27\\\\GSE175416_scRNA_seq_2022_05_27\\\\GSE175416_RAW\\\\GSM5333101_Primary_Tumor-matrix.mtx.gz"
cells_1 <- "J:\\\\Biostat_Interns\\\\Nobles_Gabrielle\\\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE175416_mouse_primary_MB_10x_2022_05_27\\GSE175416_scRNA_seq_2022_05_27\\GSE175416_RAW\\GSM5333101_Primary_Tumor-barcodes.tsv.gz"
features_1 <- "J:\\\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE175416_mouse_primary_MB_10x_2022_05_27\\GSE175416_scRNA_seq_2022_05_27\\GSE175416_RAW\\GSM5333101_Primary_Tumor-genes.tsv.gz"

# ReadMtx loads in sparse data 
counts_MB_GSE175416 <- ReadMtx(mtx_1, cells_1, features_1)


#Create a seuratobject 
seurat_object_GSE175416 = CreateSeuratObject(counts_MB_GSE175416, 
                                             project = "GSE175416")

####----GSE122871 seurat object---####

#Load in expression with Read10x
GSE122871_file<-"J:/Biostat_Interns/Nobles_Gabrielle/Moffitt_Internship_Data_2022/Internship_scRNAseq_Data_2022/Mouse_scRNAseq_Data/GSE122871_mouse_primary_GBM_10x_2022_05_27/GSE122871_scRNAseq_2022_05_27/GSE122871_barcodes_genes_matrix_2022_05_27"
list.files(GSE122871_file) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
counts_GB_GSE122871 <- Read10X(GSE122871_file) 

#Create a GSE122871 seuratobject 
seurat_object_GSE122871 <- CreateSeuratObject(counts_GB_GSE122871,
                                              project = "GSE122871")

#####----- GSE185420----#####
#Load in expression with Read10x
GSE185420_file <-"J:/Biostat_Interns/Nobles_Gabrielle/Moffitt_Internship_Data_2022/Internship_scRNAseq_Data_2022/Mouse_scRNAseq_Data/GSE185420_mouse_primary_glioblastoma_10x_Genomics_2022_05_17/GSE185420_scRNAseq_2022_05_17/GSE185420_RAW" 
list.files(GSE185420_file) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
counts_GB_GSE185420 <- Read10X(GSE185420_file)

#Create a seuratobject 
seurat_object_GSE185420 <- CreateSeuratObject(counts_GB_GSE185420, 
                                              project = "GSE185420")

#####----GSE177901----#####
# Path and file name of downloaded data - Change this to whatever your path is
GSE177901_file <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE177901_mouse_metastic_melanoma_10x_2022_05_27\\GSE177901_scRNAseq_2022_05_27\\GSE177901_scRNA-combinedEXP.txt.gz"

counts_melanoma_GSE177901 <- read.table(file = GSE177901_file, header = T, row.names=1, sep="", as.is=T)


####---Create a seuratobject----####
seurat_object_GSE177901 <- CreateSeuratObject(counts_melanoma_GSE177901,
                                              project = "GSE177901")

#####----GSE70630----####

# Path and file name of downloaded data - Change this to whatever your path is
GSE70630_file <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Human_scRNAseq_Data\\GSE70630_human_primary_OD_Smartseq2_2022_05_20\\GSE70630_scRNAseq_2022_05_20\\GSE70630_OG_processed_data_v2.txt.gz"

# Read in data as a dataframe
counts_OD_GSE70630 <- read.table(GSE70630_file)

####---Create a seurat_object_GSE70630----####
seurat_object_GSE70630 <- CreateSeuratObject(counts_OD_GSE70630, project = "GSE70630")


#####----GSE133531----####

#GSM3934450_ET_CT_12 files 
mtx_2 <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE133531_mouse_primary_10x_2022_05_24\\GSE133531_scRNAseq_2022_05_24\\GSM3934450_ET_CT_12_matrix.mtx.gz" # matrix  
cells_2 <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE133531_mouse_primary_10x_2022_05_24\\GSE133531_scRNAseq_2022_05_24\\GSM3934450_ET_CT_12_barcodes.tsv.gz" # barcode
features_2 <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE133531_mouse_primary_10x_2022_05_24\\GSE133531_scRNAseq_2022_05_24\\GSM3934450_ET_CT_12_genes.tsv.gz" # gene

# ReadMtx loads in sparse data 
counts_GSM3934450_ET_CT_12 <- ReadMtx(mtx_2, cells_2, features_2)

# seurat object GMS3934450  
seurat_object_GSM3934450 <- CreateSeuratObject(counts_GSM3934450_ET_CT_12, 
                                              project = "GMS3934450")

# GSM3934455_ET_PO_12 files 
mtx_3 <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE133531_mouse_primary_10x_2022_05_24\\GSE133531_scRNAseq_2022_05_24\\GSM3934455_ET_PO_12_matrix.mtx.gz" # matrix  
cells_3 <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE133531_mouse_primary_10x_2022_05_24\\GSE133531_scRNAseq_2022_05_24\\GSM3934455_ET_PO_12_barcodes.tsv.gz" # barcode
features_3 <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE133531_mouse_primary_10x_2022_05_24\\GSE133531_scRNAseq_2022_05_24\\GSM3934455_ET_PO_12_genes.tsv.gz" # gene
# ReadMtx loads in sparse data 
counts_GSM3934455_ET_PO_12 <- ReadMtx(mtx_3, cells_3, features_3)

# seurat object GSM3934455  
seurat_object_GSM3934455 <- CreateSeuratObject(counts_GSM3934455_ET_PO_12, 
                                               project = "GSM3934455")


# GSM3934455_ET_PO_12 files 
mtx_4 <- "J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE133531_mouse_primary_10x_2022_05_24\\GSE133531_scRNAseq_2022_05_24\\GSM3934454_C57BL6-P6-cortex_matrix.mtx.gz" # matrix  
cells_4 <-"J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE133531_mouse_primary_10x_2022_05_24\\GSE133531_scRNAseq_2022_05_24\\GSM3934454_C57BL6-P6-cortex_barcodes.tsv.gz"  # barcode
features_4 <-"J:\\Biostat_Interns\\Nobles_Gabrielle\\Moffitt_Internship_Data_2022\\Internship_scRNAseq_Data_2022\\Mouse_scRNAseq_Data\\GSE133531_mouse_primary_10x_2022_05_24\\GSE133531_scRNAseq_2022_05_24\\GSM3934454_C57BL6-P6-cortex_genes.tsv.gz"  # gene

# ReadMtx loads in sparse data 
counts_GSM3934454_C57BL6_P6_cortex <- ReadMtx(mtx_4, cells_4, features_4)

# Seurat object GSM3934454
seurat_object_GSM3934454 <- CreateSeuratObject(counts_GSM3934454_C57BL6_P6_cortex, 
                                               project = "GSM3934454")

seurat_object_GSE133531 <- merge(seurat_object_GSM3934450,
                                y = c(seurat_object_GSM3934455,seurat_object_GSM3934454),
                                add.cell.ids = c("GSM3934450","GSM3934455","GSM3934454"), 
                                project = "GSE133531")




######---- Combined_seurat_object QC----####
# Merge the seurat objects 
Combined_seurat_object <- merge(seurat_object_GSE129730,
                                y = c(seurat_object_GSE175416,seurat_object_GSE122871, 
                                      seurat_object_GSE185420, seurat_object_GSE177901,
                                      seurat_object_GSE70630, seurat_object_GSE133531),
                                add.cell.ids = c("MB_GSE129730", "MB_GSE175416", 
                                                 "GB_GSE122871","GB_GSE185420", 
                                                 "Melanoma_GSE177901", "OD_GSE70630",
                                                 "GSE133531"), 
                                project = "Combined_seurat_object")

# Replace the NA's with 0 
Combined_seurat_object@assays[["RNA"]]@data@x[is.na(Combined_seurat_object@assays[["RNA"]]@data@x)]<- 0
Combined_seurat_object@meta.data[is.na(Combined_seurat_object@meta.data)] <- 0 
Combined_seurat_object@assays[["RNA"]]@counts@x[is.na(Combined_seurat_object@assays[["RNA"]]@counts@x)] <- 0 

# Find high variable features of Combined_seurat_object
Combined_seurat_object <- FindVariableFeatures(Combined_seurat_object, selection.method = "vst", 
                                      nfeatures = 2000)

# Identify Top 10 HVG of Combined_seurat_object 
top10 <- head(VariableFeatures(Combined_seurat_object), 10)
top10

# ScaleData converts normalized gene expression to Z-score of Combined_seurat_object 
all.genes <- rownames(Combined_seurat_object)
Combined_seurat_object <- ScaleData(Combined_seurat_object, features = all.genes)


##----Combined_seurat_object Dimensional Reduction----##
# Run PCA for linear dimensional reduction 
Combined_seurat_object<- RunPCA(Combined_seurat_object,
                        features = VariableFeatures(object = Combined_seurat_object), 
                      seed.use = 42)


####----Combined_seurat_object Cluster the Cells----####
# Compute nearest neighbor grqph and SNN: shared nearest neighbor 
Combined_seurat_object <- FindNeighbors(Combined_seurat_object, dims = 1:10)

# Find clusters 
Combined_seurat_object <- FindClusters(Combined_seurat_object, resolution = 0.5)


####----Combined_seurat_object tSNE & UMAP----####
# Generate a UMAP for visualization 
Combined_seurat_object <- RunTSNE (Combined_seurat_object, dims = 1:10, check_duplicates = FALSE, 
                        seed.use = 1)

Combined_seurat_object <- RunUMAP(Combined_seurat_object, dims = 1:10, check_duplicates = FALSE,
                       seed.use = 42) 


# DimPlot use UMAP by default, Seurat clusters as identity  
DimPlot(Combined_seurat_object, reduction = "tsne",group.by = "orig.ident", seed = 1)

# Dimplot grouped by seurat clusters 
DimPlot(Combined_seurat_object,reduction = "tsne", label= T, repel = F, 
        label.size = 4, group.by = "seurat_clusters", seed = 1)



DimPlot(Combined_seurat_object, reduction = "umap",group.by = "orig.ident", seed = 1)


# 
# 
# ####----Normalized_seurat_object merge seurat object ----####
# # Normzlize GSE129730 and GSE175416 seurat objects 
# seurat_object_GSE129730_normalized <- NormalizeData(seurat_object_GSE129730,
#                                                     normalization.method = "CLR")
# 
# seurat_object_GSE175416_normalized <- NormalizeData(seurat_object_GSE175416,
#                                                     normalization.method = "CLR")
# 
# seurat_object_GSE122871_normalized <- NormalizeData(seurat_object_GSE122871,
#                                                     normalization.method = "CLR")
# 
# 
# # Merged normalized data 
# Normalized_seurat_object <- merge(seurat_object_GSE129730_normalized,
#                        y = c(seurat_object_GSE175416_normalized,seurat_object_GSE122871_normalized),
#                          add.cell.ids = c("GSE129730_normalized", "GSE175416_normalized","GSE122871_normalized"),
#                        project = "Normalized_seurat_object ",merge.data = TRUE)
# 

# Find high variable features 
# 
# Normalized_seurat_object<- FindVariableFeatures(Normalized_seurat_object, selection.method = "vst", 
#                                       nfeatures = 2000)

# Identify Top 10 HVG
top10 <- head(VariableFeatures(Normalized_seurat_object), 10)
top10

# # ScaleData converts normalized gene expression to Z-score 
# all.genes <- rownames(Normalized_seurat_object)
# Normalized_seurat_object<- ScaleData(Normalized_seurat_object, features = all.genes)
# 

# ####---- Normalized_seurat_object Dimensional Reduction----####
# # Run PCA for linear dimensional reduction 
# Normalized_seurat_object<- RunPCA(Normalized_seurat_object,
#                         features = VariableFeatures(object = Normalized_seurat_object), 
#                         seed.use = 42)
# 
# ####---- Normalized_seurat_object Cluster the Cells----####
# # Compute nearest neighbor grqph and SNN: shared nearest neighbor 
# Normalized_seurat_object <- FindNeighbors(Normalized_seurat_object, dims = 1:10)
# 
# # Find clusters 
# Normalized_seurat_object <- FindClusters(Normalized_seurat_object, resolution = 0.5)
# 

# ####----tSNE----####
# # Generate a TSNE and Umap for visualization 
# Normalized_seurat_object <- RunTSNE (Normalized_seurat_object,
#                                      dims = 1:10, check_duplicates = FALSE,
#                           seed.use = 1) 
# 
# Normalized_seurat_object <- RunUMAP(Normalized_seurat_object, 
#                                     dims = 1:10, check_duplicates = FALSE,
#                          seed.use = 42) 
# 

# # DimPlot use tsne by default, Seurat clusters as identity  
# DimPlot(Normalized_seurat_object, reduction = "tsne",
#         group.by = "orig.ident",seed = 1) 
# 
# DimPlot(Normalized_seurat_object,reduction = "tsne", label= T, repel = F, 
#         label.size = 4, group.by = "seurat_clusters", seed = 1)



#####---- Cell Annotation ----####

#Select the functions 
# hpca.ref <- celldex:: HumanPrimaryCellAtlasData()
monaco.ref <- celldex::MonacoImmuneData()
mouseRNA.ref <- celldex::MouseRNAseqData()


# Convert Seurat object to a single cell experiment for convince 
# sce.normalized <- as.SingleCellExperiment(DietSeurat(Normalized_seurat_object))


sce.combined <- as.SingleCellExperiment((DietSeurat(Combined_seurat_object)))

# #mouseRNAseq data sce.normalized 
# mouse.main <- SingleR(test = sce.normalized,assay.type.test = 1,ref = mouseRNA.ref,
#                        labels = mouseRNA.ref$label.main)
# mouse.fine <- SingleR(test = sce.normalized,assay.type.test = 1,ref = mouseRNA.ref,
#                        labels = mouseRNA.ref$label.fine)

#mouseRNAseq data sce.combined
mouse.main <- SingleR(test = sce.combined,assay.type.test = 1,ref = mouseRNA.ref,
                      labels = mouseRNA.ref$label.main)
mouse.fine <- SingleR(test = sce.combined,assay.type.test = 1,ref = mouseRNA.ref,
                      labels = mouseRNA.ref$label.fine)

# 
# # mouseRNA Normalized_seurat_object 
# Normalized_seurat_object@meta.data$mouse.main <- mouse.main$pruned.labels
# Normalized_seurat_object@meta.data$mouse.fine <- mouse.fine$pruned.labels

# mouseRNA Combined_seurat_object
Combined_seurat_object@meta.data$mouse.main <- mouse.main$pruned.labels
Combined_seurat_object@meta.data$mouse.fine <- mouse.fine$pruned.labels


# # mouseRNA Normalized_seurat_object
# Normalized_seurat_object<- SetIdent(Normalized_seurat_object, value = "mouse.fine")

# mouseRNA MB_combined
Combined_seurat_object <- SetIdent(Combined_seurat_object, value = "mouse.fine")

# ## Tsne Normalized_seurat_object
# DimPlot(Normalized_seurat_object, group.by = "mouse.fine", label = T , repel = F,
#         label.size = 1.7, reduction = "tsne" ,seed = 1) 
#   


## tsne Combined_seurat_object
DimPlot(Combined_seurat_object, group.by = "mouse.fine",label = T , repel = F,
        label.size = 1.7, reduction = "tsne", seed = 1) 




####----Differential Gene expression----####

all.markers <- FindAllMarkers(Normalized_seurat_object, only.pos = T, min.pct = 0.5, 
                              logfc.threshold = 0.5)
# To find the 
top5.markers<- as.data.frame(all.markers %>% group_by(cluster) %>% 
                                top_n(n = 5, wt = avg_log2FC))
top5.markers

# 
# Microglia.plot <- FeaturePlot(Normalized_seurat_object, 
#                               reduction = "umap", 
#                               features = c("Itgam","P2ry12", "Tmem119", "Cx3cr1"),
#                               label = FALSE, 
#                               order = TRUE,
#                               
# ) 
# 
# 
# Microglia.plot
# FeaturePlot()
# 
# 
# Astrocyte.plot <- FeaturePlot(Normalized_seurat_object, 
#                               reduction = "umap", 
#                               features = c("S100b","Slc1a3", "Aldh1l1", "Gfap"),
#                               label = FALSE, 
#                               order = TRUE,
#                               label.size = 1.60
# ) 
# 
# Astrocyte.plot
# 
# Fibroblast.plot <- FeaturePlot(Normalized_seurat_object, 
#                             reduction = "umap", 
#                             features = c("Fn1","Tnc"),
#                             label = FALSE, 
#                             order = TRUE,
#                             label.size = 1.60
# ) 
# 
# Fibroblast.plot
# 
# Tanycytes.plot <- FeaturePlot(Normalized_seurat_object, 
#                                reduction = "umap", 
#                                features = c("Gpr50","Col23a1"),
#                                label = FALSE, 
#                                order = TRUE,
#                                label.size = 1.60
#  ) 
# 
# Tanycytes.plot
# 
# Oligo.precursor.plot <- FeaturePlot(Normalized_seurat_object, 
#                               reduction = "umap", 
#                               features = c("Lhfpl3","Megf11", "Pcdh15", "Cspg4"),
#                               label = FALSE, 
#                               order = TRUE,
#                               label.size = 1.60
# ) 
# 
# Oligo.precursor.plot
# 
# oligodendrocyte.plot <-  FeaturePlot(Normalized_seurat_object, 
#                                 reduction = "umap", 
#                                 features = c("Mog","Mbp","Mag","Cldn11"),
#                                 label = FALSE, 
#                                 order = TRUE,
#                                 label.size = 1.60
# ) 
# 
# oligodendrocyte.plot
# 
# oligodendrocyte.plot2 <- FeaturePlot(Normalized_seurat_object,
#             reduction = "umap",
#             features = c("Olig1","Olig2","Olig3"),
#             label = F, 
#             order = T, 
#             label.size = 1.6)
# 
# 
# oligodendrocyte.plot2
# 
# # immunesystem.plot <-  FeaturePlot(Normalized_seurat_object, 
# #                                   reduction = "umap", 
# #                                   features = c("Il2ra","Ccr6","Cxcr3","Cd4"),
# #                                   label = FALSE, 
# #                                   order = TRUE,
# #                                   label.size = 1.60
# # ) 
# 
# # # Neural Progintor Cells 
# # npc.plot <-  FeaturePlot(Normalized_seurat_object, 
# #                                   reduction = "umap", 
# #                                   features = c("Cspg4","Rnf5","Sox2","Sox1"),
# #                                   label = FALSE, 
# #                                   order = TRUE,
# #                                   label.size = 1.60
# # ) 
# # 
# # npc.plot
# # 
# # # Neural Stem Cells 
# # nsc.plot <-  FeaturePlot(Normalized_seurat_object, 
# #                          reduction = "umap", 
# #                          features = c("Sox9","Prom1","Nes"),
# #                          label = FALSE, 
# #                          order = TRUE,
# #                          label.size = 1.60
# # ) 
# # nsc.plot
# # 
# 
# 
# # macrophage.plot <-  FeaturePlot(Normalized_seurat_object, 
# #                          reduction = "umap", 
# #                          features = c("Cd68","Cd163","Cd14"),
# #                          label = FALSE, 
# #                          order = TRUE,
# #                          label.size = 1.60
# # ) 
# # macrophage.plot
# # 
# 


# To write out table for cluster metadata 
write.table(Normalized_seurat_object@meta.data, file="C:\\Users\\80023619\\Documents\\Normalized_seurat_object.meta.data.txt", sep="\t")

# 
# 
# DimPlot(Normalized_seurat_object, reduction = "umap", group.by = "seurat_clusters") 
# 
# 
# ####### remove clusters
# # Can also provide multiple idents using the typical R syntax "c()"
# sub_obj <- subset(Normalized_seurat_object, idents = c("NPCs","aNSCs"), invert = TRUE)
# 
# 
# DimPlot(sub_obj, label = T , repel = F,
#         label.size = 2, reduction = "tsne" ,group.by = "mouse.fine") 
# 
# sub_obj2 <- subset(Combined_seurat_object, idents = c("NPCs","aNSCs"), invert = TRUE)
# 
# DimPlot(sub_obj2, label = T , repel = F,
#         label.size = 1.5, reduction = "tsne" ,group.by = "mouse.fine") 
