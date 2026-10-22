# Load libraries
library(dplyr)
library(Seurat)
library(patchwork)

myPath<-"/Users/mattia/Documents/Lezioni/SingleCell"
myProject<-"pbmc3k"
setwd(file.path(myPath,myProject))

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = file.path("filtered_gene_bc_matrices/hg19/"))

n_cells <- round(dim(pbmc.data)[2]*0.001)
n_genes <- 200 

# Initialize the Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = myProject, min.cells=n_cells, min.features=n_genes) 
pbmc

# Explore cell metadata
str(pbmc@meta.data)

# Calculate the percentage of mitochondrial counts
pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize QC metrics as a histogram
hist(pbmc$nFeature_RNA, breaks=200, xlab="nFeature_RNA", main="")
hist(pbmc$percent.mt, breaks=200, xlab="percent.mt", main="", xlim=c(0,10)) # zoom in 0-10
abline(v=5) # add vertical line

# FeatureScatter visualize feature-feature relationships 
plot1 <- FeatureScatter(pbmc, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot1 + NoLegend() 
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Apply filters
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize
pbmc <- NormalizeData(pbmc)

# Get normalized data matrix
expression_matrix<-GetAssayData(pbmc, layer="data")

