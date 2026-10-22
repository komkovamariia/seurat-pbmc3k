
# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:pca_dims)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Plot UMAP coloring by cluster
DimPlot(pbmc, group.by = "RNA_snn_res.0.5")
# highlight cluster ID
DimPlot(pbmc, group.by = "RNA_snn_res.0.5", label=TRUE, label.size=6)

# calculate number of cells in each cluster
table(pbmc$RNA_snn_res.0.5)

# try with multiple resolutions
res.try<-c(0.2,0.5,0.8)
pbmc <- FindClusters(pbmc, resolution = res.try)
DimPlot(pbmc, group.by = "RNA_snn_res.0.2", label=TRUE) + NoLegend()
DimPlot(pbmc, group.by = "RNA_snn_res.0.8", label=TRUE) + NoLegend()

################################################################################

# set a resolution
res_col<-"RNA_snn_res.0.5"

Tcells<-c("CD3D","CCR7","CD4","CD8A")
NKcells<-c("NKG7","GNLY")
Bcells<-c("CD79A","MS4A1")
Monocytes<-c("S100A9","LYZ","CD14","FCGR3A")
DC<-c("FCER1A", "CST3")
Platelet<-"PPBP"

VlnPlot(pbmc, features = Tcells, group.by=res_col, ncol=2)
VlnPlot(pbmc, features = NKcells, group.by=res_col, ncol=2)

# Dotplot of all markers
all_markers<-c(Tcells,NKcells,Bcells,Monocytes,DC,Platelet)
DotPlot(pbmc, features = all_markers, group.by=res_col) + RotatedAxis()

# change resolution
DotPlot(pbmc, features = all_markers, group.by="RNA_snn_res.0.8") + RotatedAxis()

################################################################################

# set identity classes (any column in the meta.data)
Idents(pbmc)<-res_col
# identify markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
# write results
pbmc.markers <- pbmc.markers[pbmc.markers$p_val_adj<=0.05,]
write.table(pbmc.markers,"pbmc.markers.txt",sep="\t",row.names = F)

# find all markers distinguishing cluster 1 from clusters 5
cluster1_vs_5.markers <- FindMarkers(pbmc, ident.1 = 1, ident.2 = 5)
head(cluster1_vs_5.markers, n = 5)

################################################################################

library(celldex)
library(SingleR)
library(pheatmap)

# load references
ref_BE <- BlueprintEncodeData()
ref_MI <- MonacoImmuneData() 

# extract log-normalized matrix
pbmc_matrix<-GetAssayData(pbmc)

# Run SingleR using BlueprintEncode reference and main annotation level
pred_BE <- SingleR(test=pbmc_matrix, ref=ref_BE, labels=ref_BE$label.main)
table(pred_BE$labels)

# Visualize the scores for all cells across all reference labels
plotScoreHeatmap(pred_BE)

# save labels inside metadata and plot UMAP
pbmc$BE_predictions<-pred_BE$labels
DimPlot(pbmc, group.by="BE_predictions", label=T, repel=T)

# Annotate with MI reference and main annotation level
pred_MI_main <- SingleR(test=pbmc_matrix, ref=ref_MI, labels=ref_MI$label.main)
table(pred_MI_main$labels)
pbmc$MI_predictions_main<-pred_MI_main$labels
DimPlot(pbmc, group.by="MI_predictions_main", label=T, repel=T) + NoLegend()

# Annotate with MI reference and fine annotation level
pred_MI_fine <- SingleR(test=pbmc_matrix, ref=ref_MI, labels=ref_MI$label.fine)
table(pred_MI_fine$labels)
pbmc$MI_predictions_fine<-pred_MI_fine$labels
DimPlot(pbmc, group.by="MI_predictions_fine", label=T, repel=T) + NoLegend()

# compare clusters with BE predictions
a <- table(pbmc@meta.data[,res_col], pbmc$BE_predictions)
a_freq<-a/rowSums(a)
pheatmap(t(a_freq), main="BE")

# compare clusters with MI predictions at main level
a <- table(pbmc@meta.data[,res_col], pbmc$MI_predictions_main)
a_freq<-a/rowSums(a)
pheatmap(t(a_freq), main="MI main")

# compare clusters with MI predictions at fine level
a <- table(pbmc@meta.data[,res_col], pbmc$MI_predictions_fine)
a_freq<-a/rowSums(a)
pheatmap(t(a_freq), main="MI fine")

################################################################################

# compare with QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1) + NoLegend() 

# set new identity classes
Idents(pbmc)<-res_col
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "CD4 T", "B", "CD8 T", "CD16+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# plot UMAP
DimPlot(pbmc, label = TRUE) + NoLegend()

# Calculate counts and percentage for each cell type
cell_counts <- table(Idents(pbmc))
composition <- data.frame(cluster = names(table(pbmc@meta.data[,res_col])),
                          annotation = names(cell_counts), 
                          percentage = 100 * as.numeric(cell_counts) / sum(cell_counts))
#write.table(composition,"composition.txt",sep="\t", row.names = F)

# Reproduce the default colours used by Seurat
seurat_colours <- scales::hue_pal()(nrow(composition))
names(seurat_colours) <- composition$annotation

# Bar plot
ggplot(composition, aes(x = annotation, y = percentage, fill = annotation)) +  
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = seurat_colours) + 
  labs(title = "Sample composition", x = "Cell type", y = "Cells (%)") 




