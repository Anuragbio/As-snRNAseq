library(cowplot)
library(ggplot2)
library(dplyr)
library(Seurat)
library(rlang)
library(TopKLists)
library(readxl)
library(openxlsx)
library("gplots")
#library(scAlign)
##scTransform requires considerable RAM so add this statement before proceeding to intergration.
options(future.globals.maxSize= 53687091200)
##Read data in 10x
Female_Control_R1_Data <- Read10X(data.dir="Counts_denovo/BD6CFR2/outs/filtered_feature_bc_matrix/")
Male_Control_R1_Data <- Read10X(data.dir="/Counts_denovo/BD6CMR1/outs/filtered_feature_bc_matrix/")
Female_Arsenic_R1_Data <- Read10X(data.dir="/Counts_denovo/BD6TFR1/outs/filtered_feature_bc_matrix/")
Male_Arsenic_R1_Data <- Read10X(data.dir="/Counts/BD6TMR2/outs/filtered_feature_bc_matrix/")
##Create Seurat objects
Female_Control_R1 <- CreateSeuratObject(counts=Female_Control_R1_Data,project="Female_Control",min.cells=5)

Male_Control_R1 <- CreateSeuratObject(counts=Male_Control_R1_Data,project="Male_Control",min.cells = 5)

Female_Arsenic_R1 <- CreateSeuratObject(counts=Female_Arsenic_R1_Data,project="Female_Arsenic",min.cells=5)

Male_Arsenic_R1 <- CreateSeuratObject(counts=Male_Arsenic_R1_Data,project="Male_Arsenic",min.cells = 5)

##Set identities for condition
Male_Arsenic_R1$stim <- "Arsenic"

Male_Control_R1$stim <- "Control"

Female_Arsenic_R1$stim <- "Arsenic"

Female_Control_R1$stim <- "Control"

##Set identities for gender and condition
Male_Arsenic_R1$gender_stim <- "male_Arsenic"

Male_Control_R1$gender_stim <- "male_Control"

Female_Arsenic_R1$gender_stim <- "female_Arsenic"

Female_Control_R1$gender_stim <- "female_Control"

##Set identities for sample
Male_Arsenic_R1$sample_id <- "male_Arsenic_R1"

Male_Control_R1$sample_id <- "male_Control_R1"

Female_Arsenic_R1$sample_id <- "female_Arsenic_R1"

Female_Control_R1$sample_id <- "female_Control_R1"

##Remove spurious features (low and high)

Male_Arsenic_R1 <- subset(Male_Arsenic_R1,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)

Male_Control_R1 <- subset(Male_Control_R1,subset=nFeature_RNA > 300 & nFeature_RNA <2500)

Female_Arsenic_R1 <- subset(Female_Arsenic_R1,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)

Female_Control_R1 <- subset(Female_Control_R1,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
##Normalize
Male_Arsenic_R1 <- SCTransform(Male_Arsenic_R1,verbose = FALSE,return.only.var.genes = FALSE)

Male_Control_R1 <- SCTransform(Male_Control_R1,verbose = FALSE,return.only.var.genes = FALSE)

Female_Arsenic_R1 <- SCTransform(Female_Arsenic_R1,verbose = FALSE,return.only.var.genes = FALSE)

Female_Control_R1 <- SCTransform(Female_Control_R1,verbose = FALSE,return.only.var.genes = FALSE)
##Feature selection
Brain.features <- SelectIntegrationFeatures(object.list=c(Male_Arsenic_R1,Male_Control_R1,Female_Arsenic_R1,Female_Control_R1),nfeatures=1500)
Brain.list <- PrepSCTIntegration(object.list = c(Male_Arsenic_R1,Male_Control_R1,Female_Arsenic_R1,Female_Control_R1),anchor.features = Brain.features,verbose = FALSE)


##Identify anchors
Brain.anchors <- FindIntegrationAnchors(object.list = Brain.list, normalization.method = "SCT", anchor.features = Brain.features, verbose = FALSE)
Brain.integrated <- IntegrateData(anchorset = Brain.anchors, normalization.method = "SCT", verbose = FALSE)
Brain.integrated <- RunPCA(object = Brain.integrated, verbose = FALSE)

##Break - Use elbowplot to determine how many dimensions to use for capturing as much variability as possible without using all of the axes. Look for a point of desaturation (elbow).
ElbowPlot(Brain.integrated)
##Break
##Use info from elbow plot to set number of dimensions for UMAP
Brain.integrated <- RunUMAP(Brain.integrated, reduction = "pca", dims = 1:10)
Brain.integrated <- FindNeighbors(Brain.integrated,reduction = "pca", dims = 1:10)
##Cluster cells, but first find the correct resolution by trying to apply a range of values
brain_clusters <- FindClusters(Brain.integrated,resolution = 0.4)
brain_clusters <- FindClusters(Brain.integrated,resolution = 0.5)
brain_clusters <- FindClusters(Brain.integrated,resolution = 0.6)
brain_clusters <- FindClusters(Brain.integrated,resolution = 0.7)
brain_clusters <- FindClusters(Brain.integrated,resolution = 0.8)
brain_clusters <- FindClusters(Brain.integrated,resolution = 0.9)
brain_clusters <- FindClusters(Brain.integrated,resolution = 1.0)
brain_clusters <- FindClusters(Brain.integrated,resolution = 1.1)
brain_clusters <- FindClusters(Brain.integrated,resolution = 1.2)
brain_clusters <- FindClusters(Brain.integrated,resolution = 1.3)
brain_clusters <- FindClusters(Brain.integrated,resolution = 1.4)
brain_clusters <- FindClusters(Brain.integrated,resolution = 1.5)
brain_clusters <- FindClusters(Brain.integrated,resolution = 1.6)
brain_clusters <- FindClusters(Brain.integrated,resolution = 1.7)
brain_clusters <- FindClusters(Brain.integrated,resolution = 1.8)
brain_clusters <- FindClusters(Brain.integrated,resolution = 1.9)
brain_clusters <- FindClusters(Brain.integrated,resolution = 2.0)
##Once a suitable resolution has been found, rerun the FindCluster with that value.
Brain.integrated <- FindClusters(Brain.integrated, resolution = 1.0)
Brain.integrated$Celltype <- Idents(Brain.integrated)
plots <- DimPlot(Brain.integrated, reduction= "umap", group.by = c("stim","gender_stim","Celltype","sample_id"),combine=FALSE)
DimPlot(Brain.integrated, group.by = "sample_id", pt.size = 0.001)
DimPlot(Brain.integrated, group.by = "gender_stim", pt.size = 0.001)
DimPlot(Brain.integrated, group.by = "stim", pt.size = 0.001)
DimPlot(Brain.integrated, group.by = "Celltype", pt.size = 0.001)
##Create new identities for DE analysis downstream
Brain.integrated$celltype.stim <- paste(Idents(Brain.integrated),Brain.integrated$stim,sep="_")
Brain.integrated$celltype.gender_stim <- paste(Idents(Brain.integrated),Brain.integrated$gender_stim,sep="_")
##Switch to RNA assay to find cluster markers
DefaultAssay(Brain.integrated) <- "RNA"
Brain.integrated <- NormalizeData(Brain.integrated, verbose = FALSE)
Brain.integrated_markers_all <- FindAllMarkers(Brain.integrated, min.pct=0.1, logfc.threshold = 0.25, only.pos = TRUE)
Brain.integrated_cluster_markers <- Brain.integrated_markers_all %>% group_by(cluster) %>% top_n(n=20, wt = avg_logFC) %>% print(n=20*35)
Brain.integrated_cluster_markers <- Brain.integrated_markers_all %>% group_by(cluster) #This is just a test
write.csv(Brain.integrated_cluster_markers, "/data/mackanholt_lab/achatur/SCRNA_Arsenic/Brain_integrated_cluster_markers_top20_res1_vijay_instructions.csv")
DoHeatmap(Brain.integrated, features = Brain.integrated_cluster_markers$gene) + NoLegend()
##Change labels to marker genes
new.cluster.ids <- c("0","Tm4","TmY14","Dopaminergic PAM neuron1","Dopaminergic PAM neuron2","Mip","Undetermined","Ensheathing neuropil associated glial cell","Kenyon Cell","Ensheathing glial cell","Ensheathing glial cell","Perineurial glial cell","Perineurial glia","Mi1 neurons","13","Reticular neuropil associated glial cell","T1","T1","Subperineurial Glia","Subperineurial Glia","Optic chiasma glial cell","Plasmatocyte","Kenyon cell","Tm5c","L2","Photoreceptors","Lawf2","Dm9","fat body cell","T1","Ensheathing neuropil associated glial cell","Ensheathing neuropil associated glial cell","Ensheathing neuropil associated glial cell","Ensheathing neuropil associated glial cell","Perineurial glia","Dopaminergic PAM neuron")
names(new.cluster.ids)<-levels(Brain.integrated)
Brain.integrated<-RenameIdents(Brain.integrated,new.cluster.ids)
DefaultAssay(Brain.integrated) <-"integrated"
DimPlot(Brain.integrated,label = TRUE) + NoLegend()
DimPlot(Brain.integrated,label = FALSE) + NoLegend()




##Perform DE
##Store previous biological cluster identities
Brain.integrated$cluster_ids <- Idents(Brain.integrated)
##replace with gender + condition based identities
Idents(Brain.integrated) <- "celltype.gender_stim"
##Test female first
Brain_Arsenic_DE_female_stim_C0 <- FindMarkers(Brain.integrated,ident.1 = "0_female_Arsenic",ident.2 = "0_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C1 <- FindMarkers(Brain.integrated,ident.1 = "1_female_Arsenic",ident.2 = "1_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C2 <- FindMarkers(Brain.integrated,ident.1 = "2_female_Arsenic",ident.2 = "2_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C3 <- FindMarkers(Brain.integrated,ident.1 = "3_female_Arsenic",ident.2 = "3_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C4 <- FindMarkers(Brain.integrated,ident.1 = "4_female_Arsenic",ident.2 = "4_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C5 <- FindMarkers(Brain.integrated,ident.1 = "5_female_Arsenic",ident.2 = "5_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C6 <- FindMarkers(Brain.integrated,ident.1 = "6_female_Arsenic",ident.2 = "6_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C7 <- FindMarkers(Brain.integrated,ident.1 = "7_female_Arsenic",ident.2 = "7_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C8 <- FindMarkers(Brain.integrated,ident.1 = "8_female_Arsenic",ident.2 = "8_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C9 <- FindMarkers(Brain.integrated,ident.1 = "9_female_Arsenic",ident.2 = "9_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C10 <- FindMarkers(Brain.integrated,ident.1 = "10_female_Arsenic",ident.2 = "10_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C11 <- FindMarkers(Brain.integrated,ident.1 = "11_female_Arsenic",ident.2 = "11_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C12 <- FindMarkers(Brain.integrated,ident.1 = "12_female_Arsenic",ident.2 = "12_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C13 <- FindMarkers(Brain.integrated,ident.1 = "13_female_Arsenic",ident.2 = "13_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C14 <- FindMarkers(Brain.integrated,ident.1 = "14_female_Arsenic",ident.2 = "14_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C15 <- FindMarkers(Brain.integrated,ident.1 = "15_female_Arsenic",ident.2 = "15_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C16 <- FindMarkers(Brain.integrated,ident.1 = "16_female_Arsenic",ident.2 = "16_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C17 <- FindMarkers(Brain.integrated,ident.1 = "17_female_Arsenic",ident.2 = "17_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C18 <- FindMarkers(Brain.integrated,ident.1 = "18_female_Arsenic",ident.2 = "18_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C19 <- FindMarkers(Brain.integrated,ident.1 = "19_female_Arsenic",ident.2 = "19_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C20 <- FindMarkers(Brain.integrated,ident.1 = "20_female_Arsenic",ident.2 = "20_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C21 <- FindMarkers(Brain.integrated,ident.1 = "21_female_Arsenic",ident.2 = "21_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C22 <- FindMarkers(Brain.integrated,ident.1 = "22_female_Arsenic",ident.2 = "22_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C23 <- FindMarkers(Brain.integrated,ident.1 = "23_female_Arsenic",ident.2 = "23_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C24 <- FindMarkers(Brain.integrated,ident.1 = "24_female_Arsenic",ident.2 = "24_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C25 <- FindMarkers(Brain.integrated,ident.1 = "25_female_Arsenic",ident.2 = "25_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C26 <- FindMarkers(Brain.integrated,ident.1 = "26_female_Arsenic",ident.2 = "26_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C27 <- FindMarkers(Brain.integrated,ident.1 = "27_female_Arsenic",ident.2 = "27_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C28 <- FindMarkers(Brain.integrated,ident.1 = "28_female_Arsenic",ident.2 = "28_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C29 <- FindMarkers(Brain.integrated,ident.1 = "29_female_Arsenic",ident.2 = "29_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C30 <- FindMarkers(Brain.integrated,ident.1 = "30_female_Arsenic",ident.2 = "30_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C31 <- FindMarkers(Brain.integrated,ident.1 = "31_female_Arsenic",ident.2 = "31_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C32 <- FindMarkers(Brain.integrated,ident.1 = "32_female_Arsenic",ident.2 = "32_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C33 <- FindMarkers(Brain.integrated,ident.1 = "33_female_Arsenic",ident.2 = "33_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C34 <- FindMarkers(Brain.integrated,ident.1 = "34_female_Arsenic",ident.2 = "34_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_female_stim_C35 <- FindMarkers(Brain.integrated,ident.1 = "35_female_Arsenic",ident.2 = "35_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
##test male second
Brain_Arsenic_DE_male_stim_C0 <- FindMarkers(Brain.integrated,ident.1 = "0_male_Arsenic",ident.2 = "0_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C1 <- FindMarkers(Brain.integrated,ident.1 = "1_male_Arsenic",ident.2 = "1_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C2 <- FindMarkers(Brain.integrated,ident.1 = "2_male_Arsenic",ident.2 = "2_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C3 <- FindMarkers(Brain.integrated,ident.1 = "3_male_Arsenic",ident.2 = "3_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C4 <- FindMarkers(Brain.integrated,ident.1 = "4_male_Arsenic",ident.2 = "4_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C5 <- FindMarkers(Brain.integrated,ident.1 = "5_male_Arsenic",ident.2 = "5_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C6 <- FindMarkers(Brain.integrated,ident.1 = "6_male_Arsenic",ident.2 = "6_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C7 <- FindMarkers(Brain.integrated,ident.1 = "7_male_Arsenic",ident.2 = "7_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C8 <- FindMarkers(Brain.integrated,ident.1 = "8_male_Arsenic",ident.2 = "8_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C9 <- FindMarkers(Brain.integrated,ident.1 = "9_male_Arsenic",ident.2 = "9_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C10 <- FindMarkers(Brain.integrated,ident.1 = "10_male_Arsenic",ident.2 = "10_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C11 <- FindMarkers(Brain.integrated,ident.1 = "11_male_Arsenic",ident.2 = "11_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C12 <- FindMarkers(Brain.integrated,ident.1 = "12_male_Arsenic",ident.2 = "12_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C13 <- FindMarkers(Brain.integrated,ident.1 = "13_male_Arsenic",ident.2 = "13_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C14 <- FindMarkers(Brain.integrated,ident.1 = "14_male_Arsenic",ident.2 = "14_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C15 <- FindMarkers(Brain.integrated,ident.1 = "15_male_Arsenic",ident.2 = "15_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C16 <- FindMarkers(Brain.integrated,ident.1 = "16_male_Arsenic",ident.2 = "16_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C17 <- FindMarkers(Brain.integrated,ident.1 = "17_male_Arsenic",ident.2 = "17_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C18 <- FindMarkers(Brain.integrated,ident.1 = "18_male_Arsenic",ident.2 = "18_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C19 <- FindMarkers(Brain.integrated,ident.1 = "19_male_Arsenic",ident.2 = "19_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C20 <- FindMarkers(Brain.integrated,ident.1 = "20_male_Arsenic",ident.2 = "20_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C21 <- FindMarkers(Brain.integrated,ident.1 = "21_male_Arsenic",ident.2 = "21_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C22 <- FindMarkers(Brain.integrated,ident.1 = "22_male_Arsenic",ident.2 = "22_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C23 <- FindMarkers(Brain.integrated,ident.1 = "23_male_Arsenic",ident.2 = "23_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C24 <- FindMarkers(Brain.integrated,ident.1 = "24_male_Arsenic",ident.2 = "24_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C25 <- FindMarkers(Brain.integrated,ident.1 = "25_male_Arsenic",ident.2 = "25_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C26 <- FindMarkers(Brain.integrated,ident.1 = "26_male_Arsenic",ident.2 = "26_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C27 <- FindMarkers(Brain.integrated,ident.1 = "27_male_Arsenic",ident.2 = "27_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C28 <- FindMarkers(Brain.integrated,ident.1 = "28_male_Arsenic",ident.2 = "28_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C29 <- FindMarkers(Brain.integrated,ident.1 = "29_male_Arsenic",ident.2 = "29_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C30 <- FindMarkers(Brain.integrated,ident.1 = "30_male_Arsenic",ident.2 = "30_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C31 <- FindMarkers(Brain.integrated,ident.1 = "31_male_Arsenic",ident.2 = "31_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C32 <- FindMarkers(Brain.integrated,ident.1 = "32_male_Arsenic",ident.2 = "32_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C33 <- FindMarkers(Brain.integrated,ident.1 = "33_male_Arsenic",ident.2 = "33_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C34 <- FindMarkers(Brain.integrated,ident.1 = "34_male_Arsenic",ident.2 = "34_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
Brain_Arsenic_DE_male_stim_C35 <- FindMarkers(Brain.integrated,ident.1 = "35_male_Arsenic",ident.2 = "35_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
####Saved R image until here save.image(file='/data/mackanholt_lab/achatur/SCRNA_Arsenic/SCRNA_Arsenic.RData')

# Create the list with named objects
Male_DE_list <- list(
  "Male_C0" = Brain_Arsenic_DE_male_stim_C0,
  "Male_C1" = Brain_Arsenic_DE_male_stim_C1,
  "Male_C2" = Brain_Arsenic_DE_male_stim_C2,
  "Male_C3" = Brain_Arsenic_DE_male_stim_C3,
  "Male_C4" = Brain_Arsenic_DE_male_stim_C4,
  "Male_C5" = Brain_Arsenic_DE_male_stim_C5,
  "Male_C6" = Brain_Arsenic_DE_male_stim_C6,
  "Male_C7" = Brain_Arsenic_DE_male_stim_C7,
  "Male_C8" = Brain_Arsenic_DE_male_stim_C8,
  "Male_C9" = Brain_Arsenic_DE_male_stim_C9,
  "Male_C10" = Brain_Arsenic_DE_male_stim_C10,
  "Male_C11" = Brain_Arsenic_DE_male_stim_C11,
  "Male_C12" = Brain_Arsenic_DE_male_stim_C12,
  "Male_C13" = Brain_Arsenic_DE_male_stim_C13,
  "Male_C14" = Brain_Arsenic_DE_male_stim_C14,
  "Male_C15" = Brain_Arsenic_DE_male_stim_C15,
  "Male_C16" = Brain_Arsenic_DE_male_stim_C16,
  "Male_C17" = Brain_Arsenic_DE_male_stim_C17,
  "Male_C18" = Brain_Arsenic_DE_male_stim_C18,
  "Male_C19" = Brain_Arsenic_DE_male_stim_C19,
  "Male_C20" = Brain_Arsenic_DE_male_stim_C20,
  "Male_C21" = Brain_Arsenic_DE_male_stim_C21,
  "Male_C22" = Brain_Arsenic_DE_male_stim_C22,
  "Male_C23" = Brain_Arsenic_DE_male_stim_C23,
  "Male_C24" = Brain_Arsenic_DE_male_stim_C24,
  "Male_C25" = Brain_Arsenic_DE_male_stim_C25,
  "Male_C26" = Brain_Arsenic_DE_male_stim_C26,
  "Male_C27" = Brain_Arsenic_DE_male_stim_C27,
  "Male_C28" = Brain_Arsenic_DE_male_stim_C28,
  "Male_C29" = Brain_Arsenic_DE_male_stim_C29,
  "Male_C30" = Brain_Arsenic_DE_male_stim_C30,
  "Male_C31" = Brain_Arsenic_DE_male_stim_C31
)

# Create an empty data frame to store the combined results
combined_data <- data.frame()

# Iterate through each item in the list, append a 'Cluster' column, and combine the data
for (name in names(Male_DE_list)) {
  # Get the data from the list
  data <- Male_DE_list[[name]]
  
  # Add a new column to identify the cluster
  data$Cluster <- name
  
  # Combine the data into one data frame
  combined_data <- rbind(combined_data, data)
}

# Write the combined data to a single CSV file
write.csv(combined_data, file = "/data/mackanholt_lab/achatur/SCRNA_Arsenic/Male_DE_combined_results.csv", row.names = FALSE)

cat("All results have been written to 'Male_DE_combined_results.csv'.")

