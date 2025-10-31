
library(knitr)
library(pheatmap)
library(rtracklayer)
library(chromVARmotifs)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ArchR)


set.seed(1)
addArchRThreads(threads = 1)
addArchRGenome("mm10")


load("my_session.RData")

#### ChromVar deviations

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


##motifs <- c("GATA1", "CEBPA",  "IRF4", "TBX21", "PAX5")
#*motifs <- c("Olfm4", "Lgr5", "Clu", "Muc2")
#*markerMotifs <- getFeatures(demo_proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
#*markerMotifs

#*markerMotifs <- grep("z:", markerMotifs, value = TRUE)
#*markerMotifs

#*demo_proj <- addImputeWeights(demo_proj)

#*p <- plotGroups(ArchRProj = demo_proj, 
#*  groupBy = "Clusters", 
#*  colorBy = "MotifMatrix", 
#*  name = markerMotifs,
#*  imputeWeights = getImputeWeights(demo_proj)
#*)


#*p2 <- lapply(seq_along(p), function(x){
#*  if(x != 1){
#*    p[[x]] + guides(color = FALSE, fill = FALSE) + 
#*    theme_ArchR(baseSize = 6) +
#*    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
#*    theme(
#*        axis.text.y=element_blank(), 
#*        axis.ticks.y=element_blank(),
#*        axis.title.y=element_blank()
#*    ) + ylab("")
#*  }else{
#*    p[[x]] + guides(color = FALSE, fill = FALSE) + 
#*    theme_ArchR(baseSize = 6) +
#*    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
#*    theme(
#*        axis.ticks.y=element_blank(),
#*        axis.title.y=element_blank()
#*    ) + ylab("")
#*  }
#*})


#*pdf("fastingmotif.pdf")
#*do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
#*dev.off()


#Save ArchR project
saveArchRProject(ArchRProj = demo_proj, load = TRUE)

save.image(file = "my_session12.RData")

