library(knitr)
library(pheatmap)
library(rtracklayer)
library(chromVARmotifs)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ArchR)


set.seed(1)
addArchRThreads(threads = 1)
addArchRGenome("mm10")


inputFiles <- system("ls ../*/outs/atac_fragments.tsv.gz", intern=T)
inputFiles
names(inputFiles) <- sapply(inputFiles,function(x){strsplit(x,"/")[[1]][2]})


ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles

samples <- names(inputFiles)

## Doublet identification
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1, verbose = T)

save.image(file = "fasting_atac.RData")

## Create an ArchRProject

demo_proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "results",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(demo_proj)


## Let’s plot some summary plots

pdf("fasting_fragment_size_distribution.pdf")
## plot fragment size distribution
plotFragmentSizes(ArchRProj = demo_proj)
dev.off()


## plot TSS enrichment
p1 <- plotGroups( ArchRProj = demo_proj,
 groupBy = "Sample", 
 colorBy = "cellColData",
 name = "TSSEnrichment",
 plotAs = "ridges")

pdf("fasting_peaks_TSS_enrichment.pdf")
print (p1)
dev.off()

p2 <- plotGroups( ArchRProj = demo_proj, 
groupBy = "Sample", 
colorBy = "cellColData",
name = "TSSEnrichment",
plotAs = "violin",
alpha = 0.4,
addBoxPlot = TRUE)

pdf("fasting_violin_TSS_enrichment.pdf")
print (p2)
dev.off()


pdf("fasting_TSS_enrichment_from_center.pdf")
plotTSSEnrichment(ArchRProj = demo_proj)
dev.off()



#filter out doublets
demo_proj <- filterDoublets(demo_proj)

#save the ArchR Project
saveArchRProject(ArchRProj = demo_proj,  load = FALSE)

## Dimensionality reduction & clustering

demo_proj <- addIterativeLSI(
  ArchRProj = demo_proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", iterations = 2,
  clusterParams = list( 
    resolution = 0.2, #might want to start with lower resolution. This parameter is passed to Seurat's FindClusters
    sampleCells = 5000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

demo_proj <- addClusters(
  input = demo_proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.2, #The final clustering is done in the LSI space. Once LSI is finalized, tune this parameter to control number of clusters. It is passed to Seurat's FindClusters.
  force = T
)


cMatrix <- confusionMatrix(paste0(demo_proj$Clusters), paste0(demo_proj$Sample))
cMatrix <- cMatrix / Matrix::rowSums(cMatrix)

pdf("fasting_distribution_of_clusters.pdf")
pheatmap::pheatmap(mat = cMatrix, color = paletteContinuous("blueYellow"), border_color = "black")
dev.off()

demo_proj <- addUMAP(
  ArchRProj = demo_proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP",
  nNeighbors = 30, # number of nn to compute a UMAP 
  minDist = 0.5, # how tightly UMAP can pack points
  metric = "cosine")


#Save ArchR project
saveArchRProject(ArchRProj = demo_proj, load = TRUE)

#color UMAP by “Sample”
p1 <- plotEmbedding(ArchRProj = demo_proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")


pdf("UMAP_by_sample.pdf")
print (p1)
dev.off()


#color UMAP by “Clusters”
p2 <- plotEmbedding(ArchRProj = demo_proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")


pdf("UMAP_by_cluster.pdf")
print (p2)
dev.off()


#save an editable vectorized version of this plot
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = demo_proj, addDOC = FALSE, width = 5, height = 5)    ######**************  NOTED  ***********



## Identifying Marker Genes

markersGS <- getMarkerFeatures(
  ArchRProj = demo_proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


#### We can get a list of DataFrame objects, one for each of our clusters, containing the relevant marker features using the getMarkers() function.

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList
## names(13): C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13

saveRDS(markersGS, file = paste0("results/", "fasting_markerGS.rds"))     ######**************  NOTED  ***********



#### We can visualize all of the markerList in a heatmap and label some of genes of interest. To randomly sample 40 unique gene names from the top 5 marker genes of each cluster in a markerList.
markerList.genes <- sample(unique(unlist(lapply(1:length(markerList),function(x){markerList[[names(markerList)[x]]]$name[1:5]}))),40)

heatmap.geneMK1 <- plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25",labelMarkers = markerList.genes, transpose = T)

pdf("fasting_heatmap_gene.pdf")
ComplexHeatmap::draw(heatmap.geneMK1, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


p3 <- plotEmbedding(ArchRProj = demo_proj, colorBy = "GeneScoreMatrix", name = markerList.genes, embedding = "UMAP", quantCut = c(0.01,0.95), imputeWeights = NULL)

pdf("fasting_marker_gene_UMAP_embedding.pdf")
print (p3)
dev.off()



p4 <- plotBrowserTrack(ArchRProj = demo_proj, groupBy = "Clusters", geneSymbol = markerList.genes, upstream = 50000, downstream = 50000)

plotPDF(plotList = p4, name = "clusters_marker-gene-tracks.pdf", ArchRProj = demo_proj, addDOC =F, width =5, height =5)    ########### **************** NOTED  **********

p5 <- plotBrowserTrack(ArchRProj = demo_proj, groupBy = "Sample", geneSymbol = markerList.genes, upstream = 50000, downstream = 50000)

plotPDF(plotList = p5, name = "samples_marker-gene-tracks.pdf", ArchRProj = demo_proj, addDOC =F, width =5, height =5)    ########### **************** NOTED  **********


demo_proj <- addGroupCoverages(ArchRProj = demo_proj, groupBy = "Clusters")

demo_proj <- addReproduciblePeakSet(
    ArchRProj = demo_proj,
    groupBy = "Clusters",
    pathToMacs2 = '/rsrch3/home/genomic_med/aksaw/anaconda3/envs/spatial_atac/bin/macs2', 
    #For MACS3, find this path using "which macs3" in the terminal
    #For MACS2, use pathToMacs2 = findMacs2()
    method = "q"
)

#We can similarly examine this merged peak set.
getPeakSet(demo_proj)


## Marker peaks : Often times, we are interested to know which peaks are unique to an individual cluster or a small group of clusters. These can be very useful in understanding cluster- or cell type-specific biology. We can do this in an unsupervised fashion in ArchR using the addMarkerFeatures() function in combination with useMatrix = "PeakMatrix".

demo_proj <- addPeakMatrix(demo_proj)


getAvailableMatrices(demo_proj)
## [1] "GeneScoreMatrix" "PeakMatrix"      "TileMatrix"


markersPeaks <- getMarkerFeatures(
    ArchRProj = demo_proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


markerList12 <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList12[[2]]

#visualizing marker peaks as a heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Fasting_Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = demo_proj, addDOC = FALSE)


require(BSgenome.Hsapiens.UCSC.mm10)
demo_proj <- addMotifAnnotations(ArchRProj = demo_proj, motifSet = "cisbp", name = "Motif")

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = demo_proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

#visualize using heatmap
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

plotPDF(heatmapEM, name = "Fasting_Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = demo_proj, addDOC = FALSE)

#Save ArchR project
saveArchRProject(ArchRProj = demo_proj, load = TRUE)

## ChromVar deviations

demo_proj <- addBgdPeaks(demo_proj)
## Identifying Background Peaks!
demo_proj <- addDeviationsMatrix(
  ArchRProj = demo_proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(demo_proj, name = "MotifMatrix", plot = TRUE)

pdf("fasting_variability.pdf")
plotVarDev
dev.off()


#Save ArchR project
saveArchRProject(ArchRProj = demo_proj, load = TRUE)

save.image(file = "my_session.RData")

















