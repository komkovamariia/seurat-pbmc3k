
# Calculate 2000 highly variable features
pbmc <- FindVariableFeatures(pbmc)

# Plot the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Scale data
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# Perform PCA
pbmc <- RunPCA(pbmc)
# [genes X PCs] matrix of weights
head(pbmc@reductions$pca@feature.loadings)
# [cells X PCs] matrix of cell coordinates
head(pbmc@reductions$pca@cell.embeddings)

# Examine and visualize PCA results
print(pbmc[["pca"]], nfeatures = 5)

# Visualize PCA
DimPlot(pbmc)

# Visualize features on a dimensional reduction plot
FeaturePlot(pbmc, features = "nFeature_RNA")
FeaturePlot(pbmc, features = "CD3D")

# cell cycle markers (G2/M and S phase)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# calculate scores
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes)

# visualize scores
FeaturePlot(pbmc, features = "S.Score")
FeaturePlot(pbmc, features = "G2M.Score")

# we can specify multiple features and arrange panels with ncol argument
FeaturePlot(pbmc, features = c("TOP2A","BIRC5","MKI67","PCNA"), ncol=2)

# Heatmaps of first 2 components
DimHeatmap(pbmc, dims = 1:2, cells = 500)

# Heatmaps of first 15 components saved in pdf
pdf("Heatmap_PCs.pdf", height=15)
DimHeatmap(pbmc, dims = 1:15, cells = 500)
dev.off()

# Elbow plot
ElbowPlot(pbmc)

# Specify chosen PCs number
pca_dims<-10

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:pca_dims)
str(pbmc@reductions)

# Plot UMAP embedding
DimPlot(pbmc)
DimPlot(pbmc, reduction = "pca")

# Compare cells with PPBP (platelet marker) expression in the two embeddings
FeaturePlot(pbmc, features = "PPBP")
FeaturePlot(pbmc, features = "PPBP", reduction="pca")

