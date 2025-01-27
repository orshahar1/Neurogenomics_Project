
### Genomics project part 2

################################################################
#Or Shahar 206582017
#Daniel Brooker 315015594
################################################################

library(readxl)
library(patchwork)
library(Seurat)
library(ggplot2)
library(dplyr)


############# Pre-processing


#1: Load the data
expression_data <- read_excel("C:/Users/danie/OneDrive/מסמכים/R_programing/BreastCancerExpressionInCells.xlsx")
# Transpose the matrix
expression_data <- data.frame(t(expression_data))
expression_data <- expression_data[2:nrow(expression_data),]
location_data <- read_excel("C:/Users/danie/OneDrive/מסמכים/R_programing/BreastCancerLocationOfCells.xlsx")

#2: Create Seurat object
seurat_object <- CreateSeuratObject(counts = expression_data, project = "BreastCancer")

# Add spatial information
seurat_object[["x"]] <- location_data$X
seurat_object[["y"]] <- location_data$Y

#3: Quality control
# filtering low quality cells or empty droplets
# raw data
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Filtered data 
seurat_object <- subset(seurat_object, subset = nFeature_RNA >35 & nFeature_RNA < 220  & nCount_RNA < 2500 )
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#4: Normalize data
seurat_object <- NormalizeData(seurat_object)

#5: Find variable features
# Identify the 50 most highly variable genes
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 50)
# Name the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#6: Scale data
all_genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all_genes)

#7: Run PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object),npcs=20)

# Examine and visualize PCA results a few different ways
# list of 5 
print(seurat_object[["pca"]], dims = 1:4, nfeatures = 5)
# Visualize the genes expression levels in the PCs
VizDimLoadings(seurat_object, dims = 1:4, reduction = "pca")
# Visualize the genes in the PCs dimensions
DimPlot(seurat_object, reduction = "pca") + NoLegend()
# Heat map plot of PCs
DimHeatmap(seurat_object, dims = c(1,3, 7, 9,17,20), cells = 400, balanced = TRUE)

#8: Determine the significant PCs
ElbowPlot(seurat_object)

# Create the Elbow Plot with horizontal lines 

# Extract the standard deviation of each PC (from seurat object)
pc_std_dev <- seurat_object[["pca"]]@stdev
# Create a data frame for ggplot
pc_data <- data.frame(PC = 1:length(pc_std_dev), StdDev = pc_std_dev)
ggplot(pc_data, aes(x = PC, y = StdDev)) +
  geom_point() +
  geom_line() +
  geom_segment(aes(x = 1, xend = PC, y = StdDev, yend = StdDev), linetype = "dashed", color = "blue") +  # Horizontal lines
  labs(title = "Elbow Plot with Horizontal Lines to Y-axis", x = "Principal Component", y = "Standard Deviation") +
  theme_minimal()


############# Data analyzing


#1: Cluster the cells
#KNN and louvian algorithms 
seurat_object <- FindNeighbors(seurat_object, dims = 1:12)
seurat_object <- FindClusters(seurat_object, resolution = 0.6)

#2: Run UMAP
# Visualize the clusters 
seurat_object <- RunUMAP(seurat_object, dims = 1:12)
DimPlot(seurat_object, reduction = "umap", label = TRUE) 

#3: Find markers for each cluster
# Finding differential expressed features
all_markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#"T cells" = ("CD3G", "FOXP3", "CD8A", "CD3D", "CD3E")
#"Macrophages" = ("HLA-DRA", "CD68", "CD4")
#"B cells" = ("IGHG1", "IGHG4", "IGKC", "IGHM")
#"Tumor cells" = ("EGFR", "GRB7", "ERBB2", "PGR", "CD44", "CD24", "ALDH1A3", "EPCAM", "KRT19", "KRT18", "CDH1")
#"Fibroblasts" = ("HSPG2", "SULF1")

# Shows expression probability distributions across clusters
VlnPlot(seurat_object, features = c("CD3G", "FOXP3", "CD8A", "CD3D", "CD3E"))#T cell
VlnPlot(seurat_object, features = c("HLA-DRA", "CD68", "CD4"))#Macrophages cell
VlnPlot(seurat_object, features = c("IGHG1", "IGHG4", "IGKC", "IGHM"))#B cell
VlnPlot(seurat_object, features = c("EGFR", "GRB7", "ERBB2", "PGR", "CD44", "CD24","ALDH1A3", "EPCAM", "KRT19", "KRT18", "CDH1"))#Tumor cell
VlnPlot(seurat_object, features = c("HSPG2", "SULF1"))#Fibroblasts cell

# Visualizes feature expression on a Umap plot
FeaturePlot(seurat_object, features = c("CD3G", "FOXP3", "CD8A", "CD3D", "CD3E"))#T cell
FeaturePlot(seurat_object, features = c("IGHG1", "IGHG4", "IGKC", "IGHM"))#B cell
FeaturePlot(seurat_object, features = c("EGFR", "GRB7", "ERBB2", "PGR", "CD44", "CD24","ALDH1A3", "EPCAM", "KRT19", "KRT18", "CDH1"))#Tumor cell
FeaturePlot(seurat_object, features = c("HLA-DRA", "CD68", "CD4"))#Macrophages cell
FeaturePlot(seurat_object, features = c("HSPG2", "SULF1"))#Fibroblasts cell

# List of relevant the markers for classifying the cells
specific_markers <- c("CD3G", "FOXP3", "CD8A", "CD3D", "CD3E", 
                      "HLA-DRA", "CD68", "CD4", 
                      "IGHG1", "IGHG4", "IGKC", "IGHM", 
                      "EGFR", "GRB7", "ERBB2", "PGR", "CD44", "CD24", "ALDH1A3", "EPCAM", "KRT19", "KRT18", "CDH1", 
                      "HSPG2", "SULF1")


# An additional methods to view genes expressions level in each cluster
all_markers %>%
  filter(gene %in% specific_markers) %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_specific

# DoHeatmap function generates an expression heatmap for given cells and features.
DoHeatmap(seurat_object, features = top10_specific$gene) + NoLegend()

# Finding the three most expressed genes from the specific markers in each cluster
# Highest average_log2fold_change
top3_specific_genes <- all_markers %>%
  filter(gene %in% specific_markers) %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) %>%
  arrange(cluster, -avg_log2FC)

print(top3_specific_genes)

# Classify the clusters after analyzing all the steps above
new.cluster.ids <- c("T cells","Tumor cells","B cells", "Tumor cells","Tumor cells","Macrophages","B cells" )
names(new.cluster.ids) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, new.cluster.ids)

# Plot again with the classified types
DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


############ Q1:

# Calculate percentage of immune cells
immune_cells <- WhichCells(seurat_object, idents = c("T cells", "Macrophages", "B cells"))
percent_immune <- length(immune_cells) / ncol(seurat_object) * 100

print(paste0("Percentage of immune cells: ", round(percent_immune, 2), "%"))

# Number of cells in each cell type
cell_counts <- table(Idents(seurat_object))

# Percentage of each cell type
cell_percentages <- prop.table(cell_counts) * 100

# Pie diagram 
pie(cell_percentages, 
    labels = paste0(names(cell_percentages), "\n", round(cell_percentages, 1), "%"),
    main = "Distribution of Cell Types",
    col = rainbow(length(cell_percentages)))

legend("topright", names(cell_percentages), 
       fill = rainbow(length(cell_percentages)),
       bty = "n")

print(cell_counts)
print(cell_percentages)

############ Q2:

# Define colors for immune cells and tumor cells
cell_type_colors <- c("T cells" = "blue", 
                      "B cells" = "blue", 
                      "Macrophages" = "blue",  # All immune cells are blue
                      "Tumor cells" = "red")   # Tumor cells in red

# Use the cluster type identities
seurat_object@meta.data$cell_type <- Idents(seurat_object)

# Create the spatial plot 
spatial_plot <- ggplot(seurat_object@meta.data, aes(x = x, y = y, color = cell_type)) +
  geom_point() +
  scale_color_manual(values = cell_type_colors) +  # Set colors
  theme_minimal() +
  ggtitle("Spatial distribution of cell types") +
  guides(color = guide_legend(title = "Cell Types"))

print(spatial_plot)

############ Q3:

# Analyze PD-L1 expression
seurat_object$PD_L1_expression <- GetAssayData(seurat_object, slot = "counts")["CD274",] > 0
percent_PD_L1 <- sum(seurat_object$PD_L1_expression) / ncol(seurat_object) * 100

print(paste0("Percentage of cells expressing PD-L1: ", round(percent_PD_L1, 2), "%"))

# Visualize PD-L1 expression on Umap
FeaturePlot(seurat_object, features = "CD274", pt.size = 0.5) + theme(legend.position = "right")

