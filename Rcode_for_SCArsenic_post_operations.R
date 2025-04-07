load(file='/SCRNA_Arsenic/SCRNA_Arsenic.RData')

library(cowplot)
library(ggplot2)
library(dplyr)
library(Seurat)
library(rlang)
library(TopKLists)
library(readxl)
library(openxlsx)
library("gplots")

##### Plot the histogram of cell clusters and number of cells
# Assuming Brain.integrated is your Seurat object
cluster_data <- Brain.integrated$seurat_clusters

# Create a dataframe with cluster counts
cluster_counts <- as.data.frame(table(cluster_data))
colnames(cluster_counts) <- c("Cluster", "Cell_Count")

# Plotting the histogram using ggplot2
ggplot(cluster_counts, aes(x = Cluster, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Number of Cells per Cluster",
       x = "Cluster",
       y = "Number of Cells")

#### Try to merge clusters
# First, view the cluster identities
Idents(Brain.integrated) <- "seurat_clusters"
table(Brain.integrated$seurat_clusters)
# Merge clusters 15 and 16 by assigning them a common identity
# For example, rename both as cluster 15
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "15"] <- "16"
# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)
# Find markers for cluster 15 (which now includes cells from both cluster 15 and 16)
markers_15 <- FindMarkers(Brain.integrated, ident.1 = "15", min.pct = 0.25, logfc.threshold = 0.5)
# Subset cluster 13 cells from the main Seurat object
#################################Subset cluster 13##################################################
# First, view the cluster identities
table(Brain.integrated$seurat_clusters)

Idents(Brain.integrated) <- "seurat_clusters"

# Subset cluster 13 cells from the main Seurat object
cluster_13_subset <- subset(Brain.integrated, idents = "13")
table(cluster_13_subset$seurat_clusters)
##Switch to RNA assay to find cluster markers
DefaultAssay(cluster_13_subset) <- "RNA"
# Normalize the data for the subset
cluster_13_subset <- NormalizeData(cluster_13_subset)

# Identify highly variable features
cluster_13_subset <- FindVariableFeatures(cluster_13_subset)

# Scale the data
cluster_13_subset <- ScaleData(cluster_13_subset)

# Perform PCA for dimensionality reduction
cluster_13_subset <- RunPCA(cluster_13_subset)

# Find nearest neighbors and cluster the cells
cluster_13_subset <- FindNeighbors(cluster_13_subset, dims = 1:10)  # Adjust dimensions based on PCA analysis
cluster_13_subset <- FindClusters(cluster_13_subset, resolution = 0.2)  # Adjust resolution for subclusters

# Run UMAP for visualization of subclusters
cluster_13_subset <- RunUMAP(cluster_13_subset, dims = 1:10)

# Plot UMAP to visualize subclusters
DimPlot(cluster_13_subset, label = TRUE)

cluster_13_subset$seurat_clusters <- Idents(cluster_13_subset)
plots <- DimPlot(cluster_13_subset, reduction= "umap", group.by = c("stim","gender_stim","Celltype","sample_id"),combine=FALSE)
DimPlot(cluster_13_subset, group.by = "sample_id", pt.size = 0.001)
DimPlot(cluster_13_subset, group.by = "gender_stim", pt.size = 0.001)
DimPlot(cluster_13_subset, group.by = "stim", pt.size = 0.001)
DimPlot(cluster_13_subset, group.by = "Celltype", pt.size = 0.001)

# Find markers for subclusters within cluster 13

subcluster_13_markers <- FindAllMarkers(cluster_13_subset, min.pct = 0.25, logfc.threshold = 0.25)
subcluster_13_cluster_markers <- subcluster_13_markers %>% group_by(cluster) %>% top_n(n=20, wt = avg_logFC)
# View top markers for each subcluster
table(subcluster_13_markers$cluster)


write.csv(subcluster_13_cluster_markers, "/SCRNA_Arsenic/subcluster_13_markers_top20_res0.2.csv")

##Perform DE
##Store previous biological cluster identities
cluster_13_subset$cluster_ids <- Idents(cluster_13_subset)
##replace with gender + condition based identities
Idents(cluster_13_subset) <- "celltype.gender_stim"
##Test female first
subcluster_13_0_markers <- FindAllMarkers(cluster_13_subset)
#Differential expression within subcluster 1 of cluster 13 (treated vs control)
Idents(cluster_13_subset) <- "seurat_clusters"
subcluster13_0 <- subset(cluster_13_subset, idents = "0") ###Not enough data
subcluster13_1 <- subset(cluster_13_subset, idents = "1")
subcluster13_2 <- subset(cluster_13_subset, idents = "2")
subcluster13_3 <- subset(cluster_13_subset, idents = "3")

# Perform DE analysis
Idents(subcluster13_0) <- "celltype.gender_stim"
deg_F_subcluster13_0 <- FindMarkers(subcluster13_0, ident.1 = "13_female_Arsenic", ident.2 = "13_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
deg_M_subcluster13_0 <- FindMarkers(subcluster13_0, ident.1 = "13_male_Arsenic", ident.2 = "13_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")

Idents(subcluster13_1) <- "celltype.gender_stim"
deg_F_subcluster13_1 <- FindMarkers(subcluster13_1, ident.1 = "13_female_Arsenic", ident.2 = "13_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
deg_M_subcluster13_1 <- FindMarkers(subcluster13_1, ident.1 = "13_male_Arsenic", ident.2 = "13_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")


Idents(subcluster13_2) <- "celltype.gender_stim"
deg_F_subcluster13_2 <- FindMarkers(subcluster13_2, ident.1 = "13_female_Arsenic", ident.2 = "13_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
deg_M_subcluster13_2 <- FindMarkers(subcluster13_2, ident.1 = "13_male_Arsenic", ident.2 = "13_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")

Idents(subcluster13_3) <- "celltype.gender_stim"
deg_F_subcluster13_3 <- FindMarkers(subcluster13_3, ident.1 = "13_female_Arsenic", ident.2 = "13_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
deg_M_subcluster13_3 <- FindMarkers(subcluster13_3, ident.1 = "13_male_Arsenic", ident.2 = "13_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")

# View the DEGs for subcluster 1
head(deg_M_subcluster13_0)
head(deg_M_subcluster13_1)
head(deg_M_subcluster13_2)
head(deg_M_subcluster13_3)

#################################Subset cluster 0##################################################
# First, view the cluster identities
table(Brain.integrated$seurat_clusters)

Idents(Brain.integrated) <- "seurat_clusters"

# Subset cluster 0 cells from the main Seurat object
cluster_0_subset <- subset(Brain.integrated, idents = "0")
table(cluster_0_subset$seurat_clusters)
##Switch to RNA assay to find cluster markers
DefaultAssay(cluster_0_subset) <- "RNA"
# Normalize the data for the subset
cluster_0_subset <- NormalizeData(cluster_0_subset)

# Identify highly variable features
cluster_0_subset <- FindVariableFeatures(cluster_0_subset)

# Scale the data
cluster_0_subset <- ScaleData(cluster_0_subset)

# Perform PCA for dimensionality reduction
cluster_0_subset <- RunPCA(cluster_0_subset)

# Find nearest neighbors and cluster the cells
cluster_0_subset <- FindNeighbors(cluster_0_subset, dims = 1:10)  # Adjust dimensions based on PCA analysis
cluster_0_subset <- FindClusters(cluster_0_subset, resolution = 0.2)  # Adjust resolution for subclusters

# Run UMAP for visualization of subclusters
cluster_0_subset <- RunUMAP(cluster_0_subset, dims = 1:10)

# Plot UMAP to visualize subclusters
DimPlot(cluster_0_subset, label = TRUE) +   labs(title = "Subclustering of Cluster 0")

cluster_0_subset$seurat_clusters <- Idents(cluster_0_subset)
plots <- DimPlot(cluster_0_subset, reduction= "umap", group.by = c("stim","gender_stim","Celltype","sample_id"),combine=FALSE)
DimPlot(cluster_0_subset, group.by = "sample_id", pt.size = 0.001)
DimPlot(cluster_0_subset, group.by = "gender_stim", pt.size = 0.001)
DimPlot(cluster_0_subset, group.by = "stim", pt.size = 0.001)
DimPlot(cluster_0_subset, group.by = "Celltype", pt.size = 0.001)

# Find markers for subclusters within cluster 0

subcluster_0_markers <- FindAllMarkers(cluster_0_subset, min.pct = 0.25, logfc.threshold = 0.25)
subcluster_0_cluster_markers <- subcluster_0_markers %>% group_by(cluster) %>% top_n(n=20, wt = avg_logFC)
# View top markers for each subcluster
head(subcluster_0_markers)
table(cluster_0_subset$seurat_clusters)

write.csv(subcluster_0_cluster_markers, "/SCRNA_Arsenic/subcluster_0_markers_top20_res0.2.csv")

##Perform DE
##Store previous biological cluster identities
cluster_0_subset$cluster_ids <- Idents(cluster_0_subset)
##replace with gender + condition based identities
Idents(cluster_0_subset) <- "celltype.gender_stim"
##Test female first

#Differential expression within subcluster 1 of cluster 0 (treated vs control)
Idents(cluster_0_subset) <- "seurat_clusters"
subcluster0_0 <- subset(cluster_0_subset, idents = "0")
subcluster0_1 <- subset(cluster_0_subset, idents = "1")
subcluster0_3 <- subset(cluster_0_subset, idents = "3")
subcluster0_4 <- subset(cluster_0_subset, idents = "4")
subcluster0_6 <- subset(cluster_0_subset, idents = "6")
subcluster0_8 <- subset(cluster_0_subset, idents = "8")

# Perform DE analysis
# Perform DE analysis
Idents(subcluster0_0) <- "celltype.gender_stim"
deg_F_subcluster0_0 <- FindMarkers(subcluster0_0, ident.1 = "0_female_Arsenic", ident.2 = "0_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
deg_M_subcluster0_0 <- FindMarkers(subcluster0_0, ident.1 = "0_male_Arsenic", ident.2 = "0_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")

Idents(subcluster0_1) <- "celltype.gender_stim"
deg_F_subcluster0_1 <- FindMarkers(subcluster0_1, ident.1 = "0_female_Arsenic", ident.2 = "0_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
deg_M_subcluster0_1 <- FindMarkers(subcluster0_1, ident.1 = "0_male_Arsenic", ident.2 = "0_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")

Idents(subcluster0_3) <- "celltype.gender_stim"
deg_F_subcluster0_3 <- FindMarkers(subcluster0_3, ident.1 = "0_female_Arsenic", ident.2 = "0_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
deg_M_subcluster0_3 <- FindMarkers(subcluster0_3, ident.1 = "0_male_Arsenic", ident.2 = "0_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")

Idents(subcluster0_4) <- "celltype.gender_stim"
deg_F_subcluster0_4 <- FindMarkers(subcluster0_4, ident.1 = "0_female_Arsenic", ident.2 = "0_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
deg_M_subcluster0_4 <- FindMarkers(subcluster0_4, ident.1 = "0_male_Arsenic", ident.2 = "0_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")

Idents(subcluster0_6) <- "celltype.gender_stim"
deg_F_subcluster0_6 <- FindMarkers(subcluster0_6, ident.1 = "0_female_Arsenic", ident.2 = "0_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
deg_M_subcluster0_6 <- FindMarkers(subcluster0_6, ident.1 = "0_male_Arsenic", ident.2 = "0_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")

Idents(subcluster0_8) <- "celltype.gender_stim"
deg_F_subcluster0_8 <- FindMarkers(subcluster0_8, ident.1 = "0_female_Arsenic", ident.2 = "0_female_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")
deg_M_subcluster0_8 <- FindMarkers(subcluster0_8, ident.1 = "0_male_Arsenic", ident.2 = "0_male_Control",test.use = "MAST", assay = "SCT", slot = "scale.data")

# View the DEGs for subcluster 1
head(deg_M_subcluster0_0)
head(deg_M_subcluster0_1)
head(deg_M_subcluster0_2)
head(deg_M_subcluster0_3)

write.csv(deg_F_subcluster0_0, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C0_0.csv")
write.csv(deg_M_subcluster0_0, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C0_0.csv")
write.csv(deg_F_subcluster0_1, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C0_1.csv")
write.csv(deg_M_subcluster0_1, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C0_1.csv")
write.csv(deg_F_subcluster0_3, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C0_3.csv")
write.csv(deg_M_subcluster0_3, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C0_3.csv")
write.csv(deg_F_subcluster0_4, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C0_4.csv")
write.csv(deg_M_subcluster0_4, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C0_4.csv")
write.csv(deg_F_subcluster0_6, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C0_6.csv")
write.csv(deg_M_subcluster0_6, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C0_6.csv")
write.csv(deg_F_subcluster0_8, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C0_8.csv")
write.csv(deg_M_subcluster0_8, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C0_8.csv")



###################Merge T1 neurons clusters and perform DEG analysis#################################
# Set Seurat identities to the clusters (if not already set)
Idents(Brain.integrated) <- "seurat_clusters"

# Check the current cluster distribution
table(Brain.integrated$seurat_clusters)

# Convert the cluster identities to characters to allow manipulation
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)

# Merge clusters 15, 16, and 28 by assigning them all a common identity (e.g., cluster 15)
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "16"] <- "15"
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "28"] <- "15"

# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)

# Set the Seurat identities back to the modified cluster assignments
Idents(Brain.integrated) <- Brain.integrated$seurat_clusters

# Verify that the clusters have been merged correctly
table(Idents(Brain.integrated))

# Find markers for the merged cluster 15 (which now includes cells from clusters 15, 16, and 28)
Brain.integrated_markers_15_merged_15_16_28 <- FindMarkers(Brain.integrated, ident.1 = "15", min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Brain.integrated_markers_15_merged_15_16_28, "/SCRNA_Arsenic/Brainintegrated_markers_15_merged_15_16_28.csv")


# View the top differentially expressed genes for the merged cluster 15
head(Brain.integrated_markers_15_merged_15_16_28)

Idents(Brain.integrated) <- "celltype.gender_stim"
DEG_cluster_15_male_arsenic <- FindMarkers(Brain.integrated,ident.1 = "15_male_Arsenic",ident.2 = "15_male_Control", subset.ident = "15",test.use = "MAST", assay = "SCT", slot = "scale.data")
DEG_cluster_15_female_arsenic <- FindMarkers(Brain.integrated,ident.1 = "15_female_Arsenic",ident.2 = "15_female_Control", subset.ident = "15",test.use = "MAST", assay = "SCT", slot = "scale.data")
write.csv(DEG_cluster_15_male_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_15_16_28_male_arsenic.csv")
write.csv(DEG_cluster_15_female_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_15_16_28female_arsenic.csv")


###################Merge Kenyon cell clusters and perform DEG analysis#################################
# Set Seurat identities to the clusters (if not already set)
Idents(Brain.integrated) <- "seurat_clusters"

# Check the current cluster distribution
table(Brain.integrated$seurat_clusters)

# Convert the cluster identities to characters to allow manipulation
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)

# Merge clusters 8, and 21 by assigning them all a common identity (e.g., cluster 15)
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "21"] <- "8"

# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)

# Set the Seurat identities back to the modified cluster assignments
Idents(Brain.integrated) <- Brain.integrated$seurat_clusters

# Verify that the clusters have been merged correctly
table(Idents(Brain.integrated))

# Find markers for the merged cluster 8 (which now includes cells from clusters 8 and 21)
Brain.integrated_markers_8_merged_8_21 <- FindMarkers(Brain.integrated, ident.1 = "8", min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Brain.integrated_markers_8_merged_8_21, "/SCRNA_Arsenic/Brainintegrated_markers_8_merged_8_21.csv")


# View the top differentially expressed genes for the merged cluster 8
head(Brain.integrated_markers_8_merged_8_21)

Idents(Brain.integrated) <- "celltype.gender_stim"
DEG_cluster_8_male_arsenic <- FindMarkers(Brain.integrated,ident.1 = "8_male_Arsenic",ident.2 = "8_male_Control", subset.ident = "8",test.use = "MAST", assay = "SCT", slot = "scale.data")
DEG_cluster_8_female_arsenic <- FindMarkers(Brain.integrated,ident.1 = "8_female_Arsenic",ident.2 = "8_female_Control", subset.ident = "8",test.use = "MAST", assay = "SCT", slot = "scale.data")
write.csv(DEG_cluster_8_male_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_8_21_male_arsenic.csv")
write.csv(DEG_cluster_8_female_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_8_21_female_arsenic.csv")


###################Merge Subperineurial cell clusters and perform DEG analysis#################################
# Set Seurat identities to the clusters (if not already set)
Idents(Brain.integrated) <- "seurat_clusters"

# Check the current cluster distribution
table(Brain.integrated$seurat_clusters)

# Convert the cluster identities to characters to allow manipulation
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)

# Merge clusters 17 and 18 by assigning them all a common identity (e.g., cluster 17)
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "18"] <- "17"

# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)

# Set the Seurat identities back to the modified cluster assignments
Idents(Brain.integrated) <- Brain.integrated$seurat_clusters

# Verify that the clusters have been merged correctly
table(Idents(Brain.integrated))

# Find markers for the merged cluster 17 (which now includes cells from clusters 17 and 18)
Brain.integrated_markers_17_merged_17_18 <- FindMarkers(Brain.integrated, ident.1 = "17", min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Brain.integrated_markers_17_merged_17_18, "/SCRNA_Arsenic/Brainintegrated_markers_17_merged_17_18.csv")


# View the top differentially expressed genes for the merged cluster 17
head(Brain.integrated_markers_17_merged_17_18)

Idents(Brain.integrated) <- "celltype.gender_stim"
DEG_cluster_17_male_arsenic <- FindMarkers(Brain.integrated,ident.1 = "17_male_Arsenic",ident.2 = "17_male_Control", subset.ident = "17",test.use = "MAST", assay = "SCT", slot = "scale.data")
DEG_cluster_17_female_arsenic <- FindMarkers(Brain.integrated,ident.1 = "17_female_Arsenic",ident.2 = "17_female_Control", subset.ident = "17",test.use = "MAST", assay = "SCT", slot = "scale.data")
write.csv(DEG_cluster_17_male_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_17_18_male_arsenic.csv")
write.csv(DEG_cluster_17_female_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_17_18_female_arsenic.csv")



###################Merge Perineurial cell clusters and perform DEG analysis#################################
# Set Seurat identities to the clusters (if not already set)
Idents(Brain.integrated) <- "seurat_clusters"

# Check the current cluster distribution
table(Brain.integrated$seurat_clusters)

# Convert the cluster identities to characters to allow manipulation
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)

# Merge clusters 11 and 33 by assigning them all a common identity (e.g., cluster 11)
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "33"] <- "11"

# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)

# Set the Seurat identities back to the modified cluster assignments
Idents(Brain.integrated) <- Brain.integrated$seurat_clusters

# Verify that the clusters have been merged correctly
table(Idents(Brain.integrated))

# Find markers for the merged cluster 11 (which now includes cells from clusters 11 and 33)
Brain.integrated_markers_11_merged_11_33 <- FindMarkers(Brain.integrated, ident.1 = "11", min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Brain.integrated_markers_11_merged_11_33, "/SCRNA_Arsenic/Brainintegrated_markers_11_merged_11_33.csv")


# View the top differentially expressed genes for the merged cluster 11
head(Brain.integrated_markers_11_merged_11_33)

Idents(Brain.integrated) <- "celltype.gender_stim"
DEG_cluster_11_male_arsenic <- FindMarkers(Brain.integrated,ident.1 = "11_male_Arsenic",ident.2 = "11_male_Control", subset.ident = "11",test.use = "MAST", assay = "SCT", slot = "scale.data")
DEG_cluster_11_female_arsenic <- FindMarkers(Brain.integrated,ident.1 = "11_female_Arsenic",ident.2 = "11_female_Control", subset.ident = "11",test.use = "MAST", assay = "SCT", slot = "scale.data")
write.csv(DEG_cluster_11_male_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_11_33_male_arsenic.csv")
write.csv(DEG_cluster_11_female_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_11_33_female_arsenic.csv")


###################Merge Ensheating neuropil associated glial cell clusters and perform DEG analysis#################################
# Set Seurat identities to the clusters (if not already set)
Idents(Brain.integrated) <- "seurat_clusters"

# Check the current cluster distribution
table(Brain.integrated$seurat_clusters)

# Convert the cluster identities to characters to allow manipulation
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)

# Merge clusters 7,29,30,31 and 32 by assigning them all a common identity (e.g., cluster 7)
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "29"] <- "7"
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "30"] <- "7"
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "31"] <- "7"
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "32"] <- "7"

# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)

# Set the Seurat identities back to the modified cluster assignments
Idents(Brain.integrated) <- Brain.integrated$seurat_clusters

# Verify that the clusters have been merged correctly
table(Idents(Brain.integrated))

# Find markers for the merged cluster 7 (which now includes cells from clusters 7,29,30,31 and 32)
Brain.integrated_markers_7_merged_7_29_30_31_32 <- FindMarkers(Brain.integrated, ident.1 = "7", min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Brain.integrated_markers_7_merged_7_29_30_31_32, "/SCRNA_Arsenic/Brainintegrated_markers_7_merged_7_29_30_31_32.csv")

# View the top differentially expressed genes for the merged cluster 7
head(Brain.integrated_markers_7_merged_7_29_30_31_32)

Idents(Brain.integrated) <- "celltype.gender_stim"
DEG_cluster_7_male_arsenic <- FindMarkers(Brain.integrated,ident.1 = "7_male_Arsenic",ident.2 = "7_male_Control", subset.ident = "7",test.use = "MAST", assay = "SCT", slot = "scale.data")
DEG_cluster_7_female_arsenic <- FindMarkers(Brain.integrated,ident.1 = "7_female_Arsenic",ident.2 = "7_female_Control", subset.ident = "7",test.use = "MAST", assay = "SCT", slot = "scale.data")
write.csv(DEG_cluster_7_male_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_7_29_30_31_32_male_arsenic.csv")
write.csv(DEG_cluster_7_female_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster7_29_30_31_32_female_arsenic.csv")


###################Merge Ensheating glial cell clusters and perform DEG analysis#################################
# Set Seurat identities to the clusters (if not already set)
Idents(Brain.integrated) <- "seurat_clusters"

# Check the current cluster distribution
table(Brain.integrated$seurat_clusters)

# Convert the cluster identities to characters to allow manipulation
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)

# Merge clusters 9 and 10 by assigning them all a common identity (e.g., cluster 9)
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "10"] <- "9"

# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)

# Set the Seurat identities back to the modified cluster assignments
Idents(Brain.integrated) <- Brain.integrated$seurat_clusters

# Verify that the clusters have been merged correctly
table(Idents(Brain.integrated))

# Find markers for the merged cluster 9 (which now includes cells from clusters 9 and 10)
Brain.integrated_markers_9_merged_9_10 <- FindMarkers(Brain.integrated, ident.1 = "9", min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Brain.integrated_markers_9_merged_9_10, "/SCRNA_Arsenic/Brainintegrated_markers_9_merged_9_10.csv")

# View the top differentially expressed genes for the merged cluster 9
head(Brain.integrated_markers_9_merged_9_10)

Idents(Brain.integrated) <- "celltype.gender_stim"
DEG_cluster_9_male_arsenic <- FindMarkers(Brain.integrated,ident.1 = "9_male_Arsenic",ident.2 = "9_male_Control", subset.ident = "9",test.use = "MAST", assay = "SCT", slot = "scale.data")
DEG_cluster_9_female_arsenic <- FindMarkers(Brain.integrated,ident.1 = "9_female_Arsenic",ident.2 = "9_female_Control", subset.ident = "9",test.use = "MAST", assay = "SCT", slot = "scale.data")
write.csv(DEG_cluster_9_male_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_9_10_male_arsenic.csv")
write.csv(DEG_cluster_9_female_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_9_10_female_arsenic.csv")



##Merge subclusters back to original object
DefaultAssay(Brain.integrated) <- "RNA"
DefaultAssay(cluster_13_subset) <- "RNA"
DefaultAssay(cluster_0_subset) <- "RNA"
# Ensure unique naming consistency between both objects
# Assuming cell names in both objects are unique and consistent

# Update the metadata of Brain.integrated with new subcluster identities
Brain.integrated$subcluster <- NA  # Create a new column for subcluster info
subcluster_cells <- Cells(cluster_13_subset)
subcluster_cells <- Cells(cluster_0_subset)

# Update subcluster identities in the main object based on the subcluster data
Brain.integrated$subcluster[subcluster_cells] <- cluster_13_subset$seurat_clusters
Brain.integrated$subcluster[subcluster_cells] <- cluster_0_subset$seurat_clusters### run these two cluster command one by one

# Optionally, merge cell identities to reflect both main clusters and subclusters
Brain.integrated$merged_cluster <- ifelse(is.na(Brain.integrated$subcluster),
                                          Brain.integrated$seurat_clusters,  # Retain main cluster for cells not in subcluster
                                          paste0(Brain.integrated$seurat_clusters, "_", Brain.integrated$subcluster))  # Annotate with subcluster

table(Idents(Brain.integrated))
# View the first few rows of the metadata
head(Brain.integrated@meta.data)
unique(Brain.integrated$subcluster)
unique(Brain.integrated$merged_cluster)

###################Merge Fat body cell clusters and perform DEG analysis#################################
# Set Seurat identities to the clusters (if not already set)
Idents(Brain.integrated) <- "seurat_clusters"

# Check the current cluster distribution
table(Brain.integrated$seurat_clusters)

# Convert the cluster identities to characters to allow manipulation
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)

# Merge clusters 27 and 13_3 by assigning them all a common identity (e.g., cluster 27)
Brain.integrated$seurat_clusters[Brain.integrated$merged_cluster == "13_3"] <- "27" ### Notice change in the use of merged_cluster instead of seurat_clusters and numbering also changed for subclusters of 13 from 0to 3 to 1to 4. That's why cluster 13_2 became cluster_13_3

# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)

# Set the Seurat identities back to the modified cluster assignments
Idents(Brain.integrated) <- Brain.integrated$seurat_clusters

# Verify that the clusters have been merged correctly
table(Idents(Brain.integrated))

# Find markers for the merged cluster 27 (which now includes cells from clusters 27 and 13_3)
Brain.integrated_markers_27_merged_27_13_3 <- FindMarkers(Brain.integrated, ident.1 = "27", min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Brain.integrated_markers_27_merged_27_13_3, "/SCRNA_Arsenic/Brainintegrated_markers_27_merged_27_13_3.csv")

# View the top differentially expressed genes for the merged cluster 27
head(Brain.integrated_markers_27_merged_27_13_3)

Idents(Brain.integrated) <- "celltype.gender_stim"
DEG_cluster_27_male_arsenic <- FindMarkers(Brain.integrated,ident.1 = "27_male_Arsenic",ident.2 = "27_male_Control", subset.ident = "27",test.use = "MAST", assay = "SCT", slot = "scale.data")
DEG_cluster_27_female_arsenic <- FindMarkers(Brain.integrated,ident.1 = "27_female_Arsenic",ident.2 = "27_female_Control", subset.ident = "27",test.use = "MAST", assay = "SCT", slot = "scale.data")
write.csv(DEG_cluster_27_male_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_27_13_3_male_arsenic.csv")
write.csv(DEG_cluster_27_female_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_27_13_3_female_arsenic.csv")


###################Merge optic chiasma glial clusters and perform DEG analysis#################################
# Set Seurat identities to the clusters (if not already set)
Idents(Brain.integrated) <- "seurat_clusters"

# Check the current cluster distribution
table(Brain.integrated$seurat_clusters)

# Convert the cluster identities to characters to allow manipulation
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)

# Merge clusters 19 and 13_2 by assigning them all a common identity (e.g., cluster 19)
Brain.integrated$seurat_clusters[Brain.integrated$merged_cluster == "13_2"] <- "19" ### Notice change in the use of merged_cluster instead of seurat_clusters and numbering also changed for subclusters of 13 from 0 to 3 to 1 to 4. That's why cluster 13_1 became cluster_13_2

# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)

# Set the Seurat identities back to the modified cluster assignments
Idents(Brain.integrated) <- Brain.integrated$seurat_clusters

# Verify that the clusters have been merged correctly
table(Idents(Brain.integrated))

# Find markers for the merged cluster 19 (which now includes cells from clusters 19 and 13_3)
Brain.integrated_markers_19_merged_19_13_2 <- FindMarkers(Brain.integrated, ident.1 = "19", min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Brain.integrated_markers_19_merged_19_13_2, "/SCRNA_Arsenic/Brainintegrated_markers_19_merged_19_13_2.csv")

# View the top differentially expressed genes for the merged cluster 19
head(Brain.integrated_markers_19_merged_19_13_2)

Idents(Brain.integrated) <- "celltype.gender_stim"
DEG_cluster_19_male_arsenic <- FindMarkers(Brain.integrated,ident.1 = "19_male_Arsenic",ident.2 = "19_male_Control", subset.ident = "19",test.use = "MAST", assay = "SCT", slot = "scale.data")
DEG_cluster_19_female_arsenic <- FindMarkers(Brain.integrated,ident.1 = "19_female_Arsenic",ident.2 = "19_female_Control", subset.ident = "19",test.use = "MAST", assay = "SCT", slot = "scale.data")
write.csv(DEG_cluster_19_male_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_19_13_2_male_arsenic.csv")
write.csv(DEG_cluster_19_female_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_19_13_2_female_arsenic.csv")


###################Merge photoreceptor clusters and perform DEG analysis#################################
# Set Seurat identities to the clusters (if not already set)
Idents(Brain.integrated) <- "seurat_clusters"

# Check the current cluster distribution
table(Brain.integrated$seurat_clusters)

# Convert the cluster identities to characters to allow manipulation
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)

# Merge clusters 24 and 13_4 by assigning them all a common identity (e.g., cluster 24)
Brain.integrated$seurat_clusters[Brain.integrated$merged_cluster == "13_4"] <- "24" ### Notice change in the use of merged_cluster instead of seurat_clusters and numbering also changed for subclusters of 13 from 0 to 3 to 1 to 4. That's why cluster 13_1 became cluster_13_4

# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)

# Set the Seurat identities back to the modified cluster assignments
Idents(Brain.integrated) <- Brain.integrated$seurat_clusters

# Verify that the clusters have been merged correctly
table(Idents(Brain.integrated))

# Find markers for the merged cluster 24 (which now includes cells from clusters 24 and 13_3)
Brain.integrated_markers_24_merged_24_13_4 <- FindMarkers(Brain.integrated, ident.1 = "24", min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Brain.integrated_markers_24_merged_24_13_4, "/SCRNA_Arsenic/Brainintegrated_markers_24_merged_24_13_4.csv")

# View the top differentially expressed genes for the merged cluster 24
head(Brain.integrated_markers_24_merged_24_13_4)

Idents(Brain.integrated) <- "celltype.gender_stim"
DEG_cluster_24_male_arsenic <- FindMarkers(Brain.integrated,ident.1 = "24_male_Arsenic",ident.2 = "24_male_Control", subset.ident = "24",test.use = "MAST", assay = "SCT", slot = "scale.data")
DEG_cluster_24_female_arsenic <- FindMarkers(Brain.integrated,ident.1 = "24_female_Arsenic",ident.2 = "24_female_Control", subset.ident = "24",test.use = "MAST", assay = "SCT", slot = "scale.data")
write.csv(DEG_cluster_24_male_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_24_13_4_male_arsenic.csv")
write.csv(DEG_cluster_24_female_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_24_13_4_female_arsenic.csv")

###################Merge dopaminergic PAM neuron clusters and perform DEG analysis#################################
# Set Seurat identities to the clusters (if not already set)
Idents(Brain.integrated) <- "seurat_clusters"

# Check the current cluster distribution
table(Brain.integrated$seurat_clusters)

# Convert the cluster identities to characters to allow manipulation
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)

# Merge clusters 3,4,34,0_3, and 0_6 by assigning them all a common identity (e.g., cluster 3)
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "4"] <- "3"
Brain.integrated$seurat_clusters[Brain.integrated$seurat_clusters == "34"] <- "3"
Brain.integrated$seurat_clusters[Brain.integrated$merged_cluster == "0_3"] <- "3"
Brain.integrated$seurat_clusters[Brain.integrated$merged_cluster == "0_6"] <- "3" ### Notice change in the use of merged_cluster instead of seurat_clusters and numbering also changed for subclusters of 0 from 0 to 8 to 1 to 9. That's why cluster 0_2 became cluster_0_3

# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)

# Set the Seurat identities back to the modified cluster assignments
Idents(Brain.integrated) <- Brain.integrated$seurat_clusters

# Verify that the clusters have been merged correctly
table(Idents(Brain.integrated))

# Find markers for the merged cluster 3 (which now includes cells from clusters 3,4,34,0_3 and 0_6)
Brain.integrated_markers_3_merged_3 <- FindMarkers(Brain.integrated, ident.1 = "3", min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Brain.integrated_markers_3_merged_3, "/SCRNA_Arsenic/Brainintegrated_markers_3_mergedPAMneuron.csv")

# View the top differentially expressed genes for the merged cluster 3
head(Brain.integrated_markers_3_merged_3)

Idents(Brain.integrated) <- "celltype.gender_stim"
DEG_cluster_3_male_arsenic <- FindMarkers(Brain.integrated,ident.1 = "3_male_Arsenic",ident.2 = "3_male_Control", subset.ident = "3",test.use = "MAST", assay = "SCT", slot = "scale.data")
DEG_cluster_3_female_arsenic <- FindMarkers(Brain.integrated,ident.1 = "3_female_Arsenic",ident.2 = "3_female_Control", subset.ident = "3",test.use = "MAST", assay = "SCT", slot = "scale.data")
write.csv(DEG_cluster_3_male_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_3_PAMneuron_male_arsenic.csv")
write.csv(DEG_cluster_3_female_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_3_PAMneuron_female_arsenic.csv")

###################Merge Tm5c neuron clusters and perform DEG analysis#################################
# Set Seurat identities to the clusters (if not already set)
Idents(Brain.integrated) <- "seurat_clusters"

# Check the current cluster distribution
table(Brain.integrated$seurat_clusters)

# Convert the cluster identities to characters to allow manipulation
Brain.integrated$seurat_clusters <- as.character(Brain.integrated$seurat_clusters)

# Merge clusters 22 and 0_8 by assigning them all a common identity (e.g., cluster 22)
Brain.integrated$seurat_clusters[Brain.integrated$merged_cluster == "0_8"] <- "22" ### Notice change in the use of merged_cluster instead of seurat_clusters and numbering also changed for subclusters of 0 from 0 to 8 to 1 to 9. That's why cluster 0_7 became cluster_0_8

# Convert the clusters back to a factor
Brain.integrated$seurat_clusters <- factor(Brain.integrated$seurat_clusters)

# Set the Seurat identities back to the modified cluster assignments
Idents(Brain.integrated) <- Brain.integrated$seurat_clusters

# Verify that the clusters have been merged correctly
table(Idents(Brain.integrated))

# Find markers for the merged cluster 22 (which now includes cells from clusters 22,4,224,0_22 and 0_6)
Brain.integrated_markers_22_merged_22 <- FindMarkers(Brain.integrated, ident.1 = "22", min.pct = 0.25, logfc.threshold = 0.5)
write.csv(Brain.integrated_markers_22_merged_22, "/SCRNA_Arsenic/Brainintegrated_markers_22_mergedPAMneuron.csv")

# View the top differentially expressed genes for the merged cluster 22
head(Brain.integrated_markers_22_merged_22)

Idents(Brain.integrated) <- "celltype.gender_stim"
DEG_cluster_22_male_arsenic <- FindMarkers(Brain.integrated,ident.1 = "22_male_Arsenic",ident.2 = "22_male_Control", subset.ident = "22",test.use = "MAST", assay = "SCT", slot = "scale.data")
DEG_cluster_22_female_arsenic <- FindMarkers(Brain.integrated,ident.1 = "22_female_Arsenic",ident.2 = "22_female_Control", subset.ident = "22",test.use = "MAST", assay = "SCT", slot = "scale.data")
write.csv(DEG_cluster_22_male_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_22_Tm5cneuron_male_arsenic.csv")
write.csv(DEG_cluster_22_female_arsenic, "/SCRNA_Arsenic/DEG_merged_cluster_22_Tm5cneuron_female_arsenic.csv")


#####Saving the DEG from clusters that don't need to be merged
write.csv(Brain_Arsenic_DE_male_stim_C26, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C26.csv")
write.csv(Brain_Arsenic_DE_female_stim_C26, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C26.csv")
write.csv(Brain_Arsenic_DE_male_stim_C23, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C23.csv")
write.csv(Brain_Arsenic_DE_female_stim_C23, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C23.csv")
write.csv(Brain_Arsenic_DE_male_stim_C25, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C25.csv")
write.csv(Brain_Arsenic_DE_female_stim_C25, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C25.csv")

write.csv(Brain_Arsenic_DE_male_stim_C12, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C12.csv")
write.csv(Brain_Arsenic_DE_female_stim_C12, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C12.csv")

write.csv(Brain_Arsenic_DE_male_stim_C5, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C5.csv")
write.csv(Brain_Arsenic_DE_female_stim_C5, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C5.csv")

write.csv(Brain_Arsenic_DE_male_stim_C20, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C20.csv")
write.csv(Brain_Arsenic_DE_female_stim_C20, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C20.csv")

write.csv(Brain_Arsenic_DE_male_stim_C14, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C14.csv")
write.csv(Brain_Arsenic_DE_female_stim_C14, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C14.csv")

write.csv(Brain_Arsenic_DE_male_stim_C1, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C1.csv")
write.csv(Brain_Arsenic_DE_female_stim_C1, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C1.csv")

write.csv(Brain_Arsenic_DE_male_stim_C2, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C2.csv")
write.csv(Brain_Arsenic_DE_female_stim_C2, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C2.csv")

write.csv(Brain_Arsenic_DE_male_stim_C6, "/SCRNA_Arsenic/Brain_Arsenic_DE_male_stim_C6.csv")
write.csv(Brain_Arsenic_DE_female_stim_C6, "/SCRNA_Arsenic/Brain_Arsenic_DE_female_stim_C6.csv")

save.image(file='/SCRNA_Arsenic/SCRNA_Arsenic.RData')

####Get scaled data
DefaultAssay(Brain.integrated)<-"SCT"
Idents(Brain.integrated)<-"seurat_clusters"
# list of the clusters for which data is to be fetched
clusters_to_fetch <- c(1, 2, 3, 5, 6, 7, 8, 9, 11, 12, 14, 15, 17, 19, 20, 22, 23, 24, 25, 26, 27, "0_7", "0_5", "0_2", "0_4", "0_1")

C27_cells_female<-WhichCells(Brain.integrated,idents=27)
Brain_Arsenic_DE_female_stim_27_index <- which(Brain_Arsenic_DE_female_stim_C27[,5]<=0.05)
C27_scale_data_female <- FetchData(Brain.integrated,vars = rownames(Brain_Arsenic_DE_female_stim_C27[Brain_Arsenic_DE_female_stim_27_index,]),cells=C27_cells_female,slot = "scale.data")

write.csv(C27_scale_data_female, "/SCRNA_Arsenic/C27_scale_data_female.csv")


# Set the default assay and identities for the Seurat object


# Define the clusters for which data is to be fetched for female
clusters_to_fetch <- c(1, 2, 3, 5, 6, 7, 8, 9, 11, 12, 14, 15, 17, 19, 20, 22, 23, 24, 25, 26, 27, "0_7", "0_5", "0_2", "0_4", "0_1")

# Set the default assay and identities for the Seurat object
DefaultAssay(Brain.integrated) <- "SCT"
Idents(Brain.integrated) <- "seurat_clusters"

# Fetch and save data for each cluster separately for female

# Cluster 1
C1_cells_female <- WhichCells(Brain.integrated, idents = 1)
C1_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C1[, 5] <= 0.05)
C1_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C1[C1_de_genes_index, ]), cells = C1_cells_female, slot = "scale.data")
write.csv(C1_scale_data_female, "/SCRNA_Arsenic/C1_scale_data_female.csv")

# Cluster 2
C2_cells_female <- WhichCells(Brain.integrated, idents = 2)
C2_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C2[, 5] <= 0.05)
C2_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C2[C2_de_genes_index, ]), cells = C2_cells_female, slot = "scale.data")
write.csv(C2_scale_data_female, "/SCRNA_Arsenic/C2_scale_data_female.csv")

# Cluster 3
C3_cells_female <- WhichCells(Brain.integrated, idents = 3)
C3_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C3[, 5] <= 0.05)
C3_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C3[C3_de_genes_index, ]), cells = C3_cells_female, slot = "scale.data")
write.csv(C3_scale_data_female, "/SCRNA_Arsenic/C3_scale_data_female.csv")

# Cluster 5
C5_cells_female <- WhichCells(Brain.integrated, idents = 5)
C5_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C5[, 5] <= 0.05)
C5_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C5[C5_de_genes_index, ]), cells = C5_cells_female, slot = "scale.data")
write.csv(C5_scale_data_female, "/SCRNA_Arsenic/C5_scale_data_female.csv")

# Cluster 6
C6_cells_female <- WhichCells(Brain.integrated, idents = 6)
C6_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C6[, 5] <= 0.05)
C6_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C6[C6_de_genes_index, ]), cells = C6_cells_female, slot = "scale.data")
write.csv(C6_scale_data_female, "/SCRNA_Arsenic/C6_scale_data_female.csv")

# Cluster 7
C7_cells_female <- WhichCells(Brain.integrated, idents = 7)
C7_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C7[, 5] <= 0.05)
C7_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C7[C7_de_genes_index, ]), cells = C7_cells_female, slot = "scale.data")
write.csv(C7_scale_data_female, "/SCRNA_Arsenic/C7_scale_data_female.csv")

# Cluster 8
C8_cells_female <- WhichCells(Brain.integrated, idents = 8)
C8_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C8[, 5] <= 0.05)
C8_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C8[C8_de_genes_index, ]), cells = C8_cells_female, slot = "scale.data")
write.csv(C8_scale_data_female, "/SCRNA_Arsenic/C8_scale_data_female.csv")

# Cluster 9
C9_cells_female <- WhichCells(Brain.integrated, idents = 9)
C9_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C9[, 5] <= 0.05)
C9_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C9[C9_de_genes_index, ]), cells = C9_cells_female, slot = "scale.data")
write.csv(C9_scale_data_female, "/SCRNA_Arsenic/C9_scale_data_female.csv")

# Cluster 11
C11_cells_female <- WhichCells(Brain.integrated, idents = 11)
C11_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C11[, 5] <= 0.05)
C11_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C11[C11_de_genes_index, ]), cells = C11_cells_female, slot = "scale.data")
write.csv(C11_scale_data_female, "/SCRNA_Arsenic/C11_scale_data_female.csv")

# Cluster 12
C12_cells_female <- WhichCells(Brain.integrated, idents = 12)
C12_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C12[, 5] <= 0.05)
C12_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C12[C12_de_genes_index, ]), cells = C12_cells_female, slot = "scale.data")
write.csv(C12_scale_data_female, "/SCRNA_Arsenic/C12_scale_data_female.csv")

# Cluster 14
C14_cells_female <- WhichCells(Brain.integrated, idents = 14)
C14_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C14[, 5] <= 0.05)
C14_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C14[C14_de_genes_index, ]), cells = C14_cells_female, slot = "scale.data")
write.csv(C14_scale_data_female, "/SCRNA_Arsenic/C14_scale_data_female.csv")

# Cluster 15
C15_cells_female <- WhichCells(Brain.integrated, idents = 15)
C15_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C15[, 5] <= 0.05)
C15_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C15[C15_de_genes_index, ]), cells = C15_cells_female, slot = "scale.data")
write.csv(C15_scale_data_female, "/SCRNA_Arsenic/C15_scale_data_female.csv")

# Cluster 17
C17_cells_female <- WhichCells(Brain.integrated, idents = 17)
C17_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C17[, 5] <= 0.05)
C17_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C17[C17_de_genes_index, ]), cells = C17_cells_female, slot = "scale.data")
write.csv(C17_scale_data_female, "/SCRNA_Arsenic/C17_scale_data_female.csv")

# Cluster 19
C19_cells_female <- WhichCells(Brain.integrated, idents = 19)
C19_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C19[, 5] <= 0.05)
C19_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C19[C19_de_genes_index, ]), cells = C19_cells_female, slot = "scale.data")
write.csv(C19_scale_data_female, "/SCRNA_Arsenic/C19_scale_data_female.csv")

# Cluster 20
C20_cells_female <- WhichCells(Brain.integrated, idents = 20)
C20_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C20[, 5] <= 0.05)
C20_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C20[C20_de_genes_index, ]), cells = C20_cells_female, slot = "scale.data")
write.csv(C20_scale_data_female, "/SCRNA_Arsenic/C20_scale_data_female.csv")

# Cluster 22
C22_cells_female <- WhichCells(Brain.integrated, idents = 22)
C22_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C22[, 5] <= 0.05)
C22_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C22[C22_de_genes_index, ]), cells = C22_cells_female, slot = "scale.data")
write.csv(C22_scale_data_female, "/SCRNA_Arsenic/C22_scale_data_female.csv")

# Cluster 23
C23_cells_female <- WhichCells(Brain.integrated, idents = 23)
C23_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C23[, 5] <= 0.05)
C23_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C23[C23_de_genes_index, ]), cells = C23_cells_female, slot = "scale.data")
write.csv(C23_scale_data_female, "/SCRNA_Arsenic/C23_scale_data_female.csv")

# Cluster 24
C24_cells_female <- WhichCells(Brain.integrated, idents = 24)
C24_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C24[, 5] <= 0.05)
C24_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C24[C24_de_genes_index, ]), cells = C24_cells_female, slot = "scale.data")
write.csv(C24_scale_data_female, "/SCRNA_Arsenic/C24_scale_data_female.csv")

# Cluster 25
C25_cells_female <- WhichCells(Brain.integrated, idents = 25)
C25_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C25[, 5] <= 0.05)
C25_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C25[C25_de_genes_index, ]), cells = C25_cells_female, slot = "scale.data")
write.csv(C25_scale_data_female, "/SCRNA_Arsenic/C25_scale_data_female.csv")

# Cluster 26
C26_cells_female <- WhichCells(Brain.integrated, idents = 26)
C26_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C26[, 5] <= 0.05)
C26_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C26[C26_de_genes_index, ]), cells = C26_cells_female, slot = "scale.data")
write.csv(C26_scale_data_female, "/SCRNA_Arsenic/C26_scale_data_female.csv")

# Cluster 27
C27_cells_female <- WhichCells(Brain.integrated, idents = 27)
C27_de_genes_index <- which(Brain_Arsenic_DE_female_stim_C27[, 5] <= 0.05)
C27_scale_data_female <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_female_stim_C27[C27_de_genes_index, ]), cells = C27_cells_female, slot = "scale.data")
write.csv(C27_scale_data_female, "/SCRNA_Arsenic/C27_scale_data_female.csv")

Idents(Brain.integrated) <- "merged_cluster" #Remember for cluster 0, 0_0 became 0_1 and so on
# Cluster 0_7
C0_7_cells_female <- WhichCells(Brain.integrated, idents = "0_7")
C0_7_de_genes_index <- which(deg_F_subcluster0_6[, 5] <= 0.05)
C0_7_scale_data_female <- FetchData(Brain.integrated, vars = rownames(deg_F_subcluster0_6[C0_7_de_genes_index, ]), cells = C0_7_cells_female, slot = "scale.data")
write.csv(C0_7_scale_data_female, "/SCRNA_Arsenic/C0_7_scale_data_female.csv")

# Cluster 0_5
C0_5_cells_female <- WhichCells(Brain.integrated, idents = "0_5")
C0_5_de_genes_index <- which(deg_F_subcluster0_4[, 5] <= 0.05)
C0_5_scale_data_female <- FetchData(Brain.integrated, vars = rownames(deg_F_subcluster0_4[C0_5_de_genes_index, ]), cells = C0_5_cells_female, slot = "scale.data")
write.csv(C0_5_scale_data_female, "/SCRNA_Arsenic/C0_5_scale_data_female.csv")

# Cluster 0_2
C0_2_cells_female <- WhichCells(Brain.integrated, idents = "0_2")
C0_2_de_genes_index <- which(deg_F_subcluster0_1[, 5] <= 0.05)
C0_2_scale_data_female <- FetchData(Brain.integrated, vars = rownames(deg_F_subcluster0_1[C0_2_de_genes_index, ]), cells = C0_2_cells_female, slot = "scale.data")
write.csv(C0_2_scale_data_female, "/SCRNA_Arsenic/C0_2_scale_data_female.csv")

# Cluster 0_4
C0_4_cells_female <- WhichCells(Brain.integrated, idents = "0_4")
C0_4_de_genes_index <- which(deg_F_subcluster0_3[, 5] <= 0.05)
C0_4_scale_data_female <- FetchData(Brain.integrated, vars = rownames(deg_F_subcluster0_3[C0_4_de_genes_index, ]), cells = C0_4_cells_female, slot = "scale.data")
write.csv(C0_4_scale_data_female, "/SCRNA_Arsenic/C0_4_scale_data_female.csv")

# Cluster 0_1
C0_1_cells_female <- WhichCells(Brain.integrated, idents = "0_1")
C0_1_de_genes_index <- which(deg_F_subcluster0_0[, 5] <= 0.05)
C0_1_scale_data_female <- FetchData(Brain.integrated, vars = rownames(deg_F_subcluster0_0[C0_1_de_genes_index, ]), cells = C0_1_cells_female, slot = "scale.data")
write.csv(C0_1_scale_data_female, "/SCRNA_Arsenic/C0_1_scale_data_female.csv")


# Define the clusters for which data is to be fetched for male
clusters_to_fetch <- c(1, 2, 3, 5, 6, 7, 8, 9, 11, 12, 14, 15, 17, 19, 20, 22, 23, 24, 25, 26, 27, "0_6", "0_4", "0_1", "0_3", "0_0")

# Set the default assay and identities for the Seurat object
DefaultAssay(Brain.integrated) <- "SCT"
Idents(Brain.integrated) <- "seurat_clusters"

# Fetch and save data for each cluster separately for male

# Cluster 1
C1_cells_male <- WhichCells(Brain.integrated, idents = 1)
C1_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C1[, 5] <= 0.05)
C1_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C1[C1_de_genes_index, ]), cells = C1_cells_male, slot = "scale.data")
write.csv(C1_scale_data_male, "/SCRNA_Arsenic/C1_scale_data_male.csv")

# Cluster 2
C2_cells_male <- WhichCells(Brain.integrated, idents = 2)
C2_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C2[, 5] <= 0.05)
C2_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C2[C2_de_genes_index, ]), cells = C2_cells_male, slot = "scale.data")
write.csv(C2_scale_data_male, "/SCRNA_Arsenic/C2_scale_data_male.csv")

# Cluster 3
C3_cells_male <- WhichCells(Brain.integrated, idents = 3)
C3_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C3[, 5] <= 0.05)
C3_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C3[C3_de_genes_index, ]), cells = C3_cells_male, slot = "scale.data")
write.csv(C3_scale_data_male, "/SCRNA_Arsenic/C3_scale_data_male.csv")

# Cluster 5
C5_cells_male <- WhichCells(Brain.integrated, idents = 5)
C5_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C5[, 5] <= 0.05)
C5_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C5[C5_de_genes_index, ]), cells = C5_cells_male, slot = "scale.data")
write.csv(C5_scale_data_male, "/SCRNA_Arsenic/C5_scale_data_male.csv")

# Cluster 6
C6_cells_male <- WhichCells(Brain.integrated, idents = 6)
C6_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C6[, 5] <= 0.05)
C6_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C6[C6_de_genes_index, ]), cells = C6_cells_male, slot = "scale.data")
write.csv(C6_scale_data_male, "/SCRNA_Arsenic/C6_scale_data_male.csv")

# Cluster 7
C7_cells_male <- WhichCells(Brain.integrated, idents = 7)
C7_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C7[, 5] <= 0.05)
C7_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C7[C7_de_genes_index, ]), cells = C7_cells_male, slot = "scale.data")
write.csv(C7_scale_data_male, "/SCRNA_Arsenic/C7_scale_data_male.csv")

# Cluster 8
C8_cells_male <- WhichCells(Brain.integrated, idents = 8)
C8_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C8[, 5] <= 0.05)
C8_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C8[C8_de_genes_index, ]), cells = C8_cells_male, slot = "scale.data")
write.csv(C8_scale_data_male, "/SCRNA_Arsenic/C8_scale_data_male.csv")

# Cluster 9
C9_cells_male <- WhichCells(Brain.integrated, idents = 9)
C9_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C9[, 5] <= 0.05)
C9_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C9[C9_de_genes_index, ]), cells = C9_cells_male, slot = "scale.data")
write.csv(C9_scale_data_male, "/SCRNA_Arsenic/C9_scale_data_male.csv")

# Cluster 11
C11_cells_male <- WhichCells(Brain.integrated, idents = 11)
C11_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C11[, 5] <= 0.05)
C11_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C11[C11_de_genes_index, ]), cells = C11_cells_male, slot = "scale.data")
write.csv(C11_scale_data_male, "/SCRNA_Arsenic/C11_scale_data_male.csv")

# Cluster 12
C12_cells_male <- WhichCells(Brain.integrated, idents = 12)
C12_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C12[, 5] <= 0.05)
C12_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C12[C12_de_genes_index, ]), cells = C12_cells_male, slot = "scale.data")
write.csv(C12_scale_data_male, "/SCRNA_Arsenic/C12_scale_data_male.csv")

# Cluster 14
C14_cells_male <- WhichCells(Brain.integrated, idents = 14)
C14_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C14[, 5] <= 0.05)
C14_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C14[C14_de_genes_index, ]), cells = C14_cells_male, slot = "scale.data")
write.csv(C14_scale_data_male, "/SCRNA_Arsenic/C14_scale_data_male.csv")

# Cluster 15
C15_cells_male <- WhichCells(Brain.integrated, idents = 15)
C15_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C15[, 5] <= 0.05)
C15_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C15[C15_de_genes_index, ]), cells = C15_cells_male, slot = "scale.data")
write.csv(C15_scale_data_male, "/SCRNA_Arsenic/C15_scale_data_male.csv")

# Cluster 17
C17_cells_male <- WhichCells(Brain.integrated, idents = 17)
C17_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C17[, 5] <= 0.05)
C17_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C17[C17_de_genes_index, ]), cells = C17_cells_male, slot = "scale.data")
write.csv(C17_scale_data_male, "/SCRNA_Arsenic/C17_scale_data_male.csv")

# Cluster 19
C19_cells_male <- WhichCells(Brain.integrated, idents = 19)
C19_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C19[, 5] <= 0.05)
C19_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C19[C19_de_genes_index, ]), cells = C19_cells_male, slot = "scale.data")
write.csv(C19_scale_data_male, "/SCRNA_Arsenic/C19_scale_data_male.csv")

# Cluster 20
C20_cells_male <- WhichCells(Brain.integrated, idents = 20)
C20_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C20[, 5] <= 0.05)
C20_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C20[C20_de_genes_index, ]), cells = C20_cells_male, slot = "scale.data")
write.csv(C20_scale_data_male, "/SCRNA_Arsenic/C20_scale_data_male.csv")

# Cluster 22
C22_cells_male <- WhichCells(Brain.integrated, idents = 22)
C22_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C22[, 5] <= 0.05)
C22_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C22[C22_de_genes_index, ]), cells = C22_cells_male, slot = "scale.data")
write.csv(C22_scale_data_male, "/SCRNA_Arsenic/C22_scale_data_male.csv")

# Cluster 23
C23_cells_male <- WhichCells(Brain.integrated, idents = 23)
C23_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C23[, 5] <= 0.05)
C23_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C23[C23_de_genes_index, ]), cells = C23_cells_male, slot = "scale.data")
write.csv(C23_scale_data_male, "/SCRNA_Arsenic/C23_scale_data_male.csv")

# Cluster 24
C24_cells_male <- WhichCells(Brain.integrated, idents = 24)
C24_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C24[, 5] <= 0.05)
C24_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C24[C24_de_genes_index, ]), cells = C24_cells_male, slot = "scale.data")
write.csv(C24_scale_data_male, "/SCRNA_Arsenic/C24_scale_data_male.csv")

# Cluster 25
C25_cells_male <- WhichCells(Brain.integrated, idents = 25)
C25_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C25[, 5] <= 0.05)
C25_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C25[C25_de_genes_index, ]), cells = C25_cells_male, slot = "scale.data")
write.csv(C25_scale_data_male, "/SCRNA_Arsenic/C25_scale_data_male.csv")

# Cluster 26
C26_cells_male <- WhichCells(Brain.integrated, idents = 26)
C26_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C26[, 5] <= 0.05)
C26_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C26[C26_de_genes_index, ]), cells = C26_cells_male, slot = "scale.data")
write.csv(C26_scale_data_male, "/SCRNA_Arsenic/C26_scale_data_male.csv")

# Cluster 27
C27_cells_male <- WhichCells(Brain.integrated, idents = 27)
C27_de_genes_index <- which(Brain_Arsenic_DE_male_stim_C27[, 5] <= 0.05)
C27_scale_data_male <- FetchData(Brain.integrated, vars = rownames(Brain_Arsenic_DE_male_stim_C27[C27_de_genes_index, ]), cells = C27_cells_male, slot = "scale.data")
write.csv(C27_scale_data_male, "/SCRNA_Arsenic/C27_scale_data_male.csv")


Idents(Brain.integrated) <- "merged_cluster" #Remember for cluster 0, 0_0 became 0_1 and so on
# Cluster 0_7
C0_7_cells_male <- WhichCells(Brain.integrated, idents = "0_7")
C0_7_de_genes_index <- which(deg_F_subcluster0_6[, 5] <= 0.05)
C0_7_scale_data_male <- FetchData(Brain.integrated, vars = rownames(deg_F_subcluster0_6[C0_7_de_genes_index, ]), cells = C0_7_cells_male, slot = "scale.data")
write.csv(C0_7_scale_data_male, "/SCRNA_Arsenic/C0_7_scale_data_male.csv")

# Cluster 0_5
C0_5_cells_male <- WhichCells(Brain.integrated, idents = "0_5")
C0_5_de_genes_index <- which(deg_F_subcluster0_4[, 5] <= 0.05)
C0_5_scale_data_male <- FetchData(Brain.integrated, vars = rownames(deg_F_subcluster0_4[C0_5_de_genes_index, ]), cells = C0_5_cells_male, slot = "scale.data")
write.csv(C0_5_scale_data_male, "/SCRNA_Arsenic/C0_5_scale_data_male.csv")

# Cluster 0_2
C0_2_cells_male <- WhichCells(Brain.integrated, idents = "0_2")
C0_2_de_genes_index <- which(deg_F_subcluster0_1[, 5] <= 0.05)
C0_2_scale_data_male <- FetchData(Brain.integrated, vars = rownames(deg_F_subcluster0_1[C0_2_de_genes_index, ]), cells = C0_2_cells_male, slot = "scale.data")
write.csv(C0_2_scale_data_male, "/SCRNA_Arsenic/C0_2_scale_data_male.csv")

# Cluster 0_4
C0_4_cells_male <- WhichCells(Brain.integrated, idents = "0_4")
C0_4_de_genes_index <- which(deg_F_subcluster0_3[, 5] <= 0.05)
C0_4_scale_data_male <- FetchData(Brain.integrated, vars = rownames(deg_F_subcluster0_3[C0_4_de_genes_index, ]), cells = C0_4_cells_male, slot = "scale.data")
write.csv(C0_4_scale_data_male, "/SCRNA_Arsenic/C0_4_scale_data_male.csv")

# Cluster 0_1
C0_1_cells_male <- WhichCells(Brain.integrated, idents = "0_1")
C0_1_de_genes_index <- which(deg_F_subcluster0_0[, 5] <= 0.05)
C0_1_scale_data_male <- FetchData(Brain.integrated, vars = rownames(deg_F_subcluster0_0[C0_1_de_genes_index, ]), cells = C0_1_cells_male, slot = "scale.data")
write.csv(C0_1_scale_data_male, "/SCRNA_Arsenic/C0_1_scale_data_male.csv")

save.image(file='/SCRNA_Arsenic/SCRNA_Arsenic.RData')




