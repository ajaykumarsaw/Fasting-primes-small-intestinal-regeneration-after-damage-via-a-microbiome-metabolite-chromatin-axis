library(Seurat)
library(ggplot2)
library(tidyverse)
#library(gridExtra)
library(patchwork)
set.seed(4567)

cat("Step 1 ##############") 

dirs <- Sys.glob("filtered_feature_bc_matrix*")
dirs
                      
for(x in dirs){                                 
  name <- gsub('filtered_feature_bc_matrix_','', x)
  
  cts <- ReadMtx(mtx = paste0(getwd(),'/',x,'/matrix.mtx.gz'),
          features = paste0(getwd(),'/',x,'/features.tsv.gz'),
          cells = paste0(getwd(),'/',x,'/barcodes.tsv.gz'))

          # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}

ls()

ls()[3:10]

cat("Step 2 ##############") 

merged_seurat <- merge(Fast_0h, y = c(Fast_24hpi, Fast_48hpi, Fast_72hpi, Fed_0h,  Fed_24hpi,  Fed_48hpi,  Fed_72hpi), add.cell.ids = c("Fast_0h", "Fast_24hpi", "Fast_48hpi", "Fast_72hpi", "Fed_0h", "Fed_24hpi", "Fed_48hpi", "Fed_72hpi"), project = "fasting_rna")

#### NOTE: sequential order must match with "add.cell.id" sequential order for correct annotation 
## add layers

merged_seurat <- JoinLayers(object = merged_seurat)

merged_seurat

cat("Step 3 ##############") 
# QC & filtering -----------------------

head(merged_seurat@meta.data)

cat("Step 4 ##############") 
# create a sample column

merged_seurat$sample <- rownames(merged_seurat@meta.data)

cat("Step 5 ##############") 
# split sample column

merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('type','time', 'Barcode'), sep = '_')   ##############**************


merged_seurat@meta.data$ID <- paste(merged_seurat@meta.data$type, merged_seurat@meta.data$time, sep="_")

head(merged_seurat@meta.data)

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')


unique(merged_seurat@meta.data$ID)

# explore QC plot IDwise indivutally

cat("Step 6 ##############") 

nab <- list("Fast_0h", "Fast_24hpi", "Fast_48hpi", "Fast_72hpi", "Fed_0h", "Fed_24hpi", "Fed_48hpi", "Fed_72hpi")        ######*************************

for(r1 in nab){
qc_1 <- VlnPlot(subset(x = merged_seurat, subset = ID == r1), features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)

qc_2 <- FeatureScatter(subset(x = merged_seurat, subset = ID == r1), feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

pdf(paste0("volcano_plot_",r1, ".pdf"))
print(qc_1)
dev.off()

pdf(paste0("scatter_plot_",r1, ".pdf"))
print(qc_2)
dev.off()
}

qc_1 <- VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
qc_11 <- VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3,  pt.size = 0)

qc_2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

pdf("volcano_plot_all_merged_samples.pdf")
print(qc_1)
dev.off()


pdf("volcano_plot_all_merged_samples_no_dot.pdf")
print(qc_11)
dev.off()


pdf("scatter_plot_all_merged_samples.pdf")
print(qc_2)
dev.off()

cat("Step 7 ##############") 


# filtering
## merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 200 & nFeature_RNA > 400 & mitoPercent < 20)

merged_seurat_filtered  <- merged_seurat

merged_seurat


qc_fil1 <- VlnPlot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
qc_fil11 <- VlnPlot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3,  pt.size = 0)

qc_fil2 <- FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

pdf("filtered_volcano_plot_all_merged_samples.pdf")
print(qc_fil1)
dev.off()


pdf("filtered_volcano_plot_all_merged_samples_no_dot.pdf")
print(qc_fil11)
dev.off()


pdf("filtered_scatter_plot_all_merged_samples.pdf")
print(qc_fil2)
dev.off()



# perform standard workflow steps to figure out if we see any batch effects --------

merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)


pdf("merged_samples_all_elbow.pdf")
ElbowPlot(merged_seurat_filtered) 
dev.off()

merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)



p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'ID')
# p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type')


#pdf("beforegrid.pdf")
#grid.arrange(p1, p2, ncol = 2, nrow = 2) 
#dev.off()

pdf("merged_samples_all_umap.pdf")
print(p1)
dev.off()


# split the dataset into a list of six seurat objects ()
obj.list <- SplitObject(merged_seurat_filtered, split.by = "ID")

# normalize and identify variable features for each dataset independently
obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})


###  saveRDS(obj.list, file = "obj_list.rds")
###  obj.list <- readRDS("obj_list.rds")


# select features that are repeatedly variable across datasets for integration

features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

### integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seurat.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering

# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)
seurat.integrated <- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:50)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.2)


p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'ID')
#p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Type')
p5 <- DimPlot(seurat.integrated, reduction = "umap", label = TRUE, repel = TRUE)


pdf("after_batcheffect_correction_merged_samples_all_umap_patientID.pdf")
print(p3)
dev.off()

pdf("after_batcheffect_correction_merged_samples_all_umap_cell_type_clustering.pdf")
print(p5)
dev.off()

pdf("after_batcheffect_correction_merged_samples_all_umap_patientID_cell_type_clustering.pdf")
p3 + p5
dev.off()


p6 <- DimPlot(seurat.integrated, reduction = "umap", split.by = "ID",label = TRUE, repel = TRUE)

pdf("after_batcheffect_correction_merged_samples_all_umap_cell_type_against_each_ID.pdf")
p6
dev.off()

# For performing differential expression after integration, we switch back to the original
# data



write.table(Idents(seurat.integrated),"Cell_with_clustering.txt",sep="\t",row.names=TRUE)

saveRDS(seurat.integrated, file = "seurat_integrated.rds")
### seurat.integrated <- readRDS("seurat_integrated.rds")

head(Idents(seurat.integrated))  ## for number of clustering identification

save.image(file="fasting.RData")


