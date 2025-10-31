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


