# https://bioconductor.org/books/3.12/OSCA/hca-human-bone-marrow-10x-genomics.html

# Setup path
figures <- "/mnt/d12b/SINGLECELL/figures/"
data <- "/mnt/d12b/SINGLECELL/testdata/"

library(tidyverse)

# Load and prepare raw sce data
## Load data from the Human Cell Atlas
library(HCAData)
sce.bone <- HCAData('ica_bone_marrow')
sce.bone$Donor <- sub("_.*", "", sce.bone$Barcode)
## We use symbols in place of IDs for easier interpretation later
library(EnsDb.Hsapiens.v86)
## Map chromosome id
rowData(sce.bone)$Chr <- mapIds(EnsDb.Hsapiens.v86, keys=rownames(sce.bone),
    column="SEQNAME", keytype="GENEID")
## Map the rownames with feature name from EnsDb
library(scater)
rownames(sce.bone) <- uniquifyFeatureNames(rowData(sce.bone)$ID,
    names = rowData(sce.bone)$Symbol)

# Write the raw data
saveRDS(sce.bone, file = paste0(data, "sce_bone_raw.rds"))

#| ----
# Quality control
## Multi-core processing
library(BiocParallel)
bpp <- BiocParallel::MulticoreParam(24)
## Add quality metrics
sce.bone <- unfiltered <- addPerCellQC(sce.bone, BPPARAM=bpp,
    subsets=list(Mito=which(rowData(sce.bone)$Chr=="MT")))

## quick quality check
qc <- quickPerCellQC(colData(sce.bone), batch=sce.bone$Donor,
    sub.fields="subsets_Mito_percent")
df <- qc |> as.data.frame() |> colSums()
saveRDS(df, file = paste0(data, "qc_metrics.rds"))

## Filter cell
sce.bone <- sce.bone[,!qc$discard]
unfiltered$discard <- qc$discard

## save data
saveRDS(unfiltered, file=paste0(data, "sce_bone_unfiltered.rds"))
saveRDS(sce.bone, file=paste0(data, "sce_bone_filtered.rds"))


# Metadata (colData)
df <- colData(sce.bone) |> as.data.frame()
### Update column name
colnames(df)
names(df)[names(df) == "sum"] <- "nUMI"
names(df)[names(df) == "detected"] <- "nGene"
### save metadata
saveRDS(df, file=paste0(data, "sce_bone_metadata.rds"))


# Plot quality data
p <- gridExtra::grid.arrange(
    plotColData(unfiltered, x="Donor", y="sum", colour_by="discard") +
        scale_y_log10() + ggtitle("Total count") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    plotColData(unfiltered, x="Donor", y="detected", colour_by="discard") +
        scale_y_log10() + ggtitle("Detected features") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    plotColData(unfiltered, x="Donor", y="subsets_Mito_percent",
        colour_by="discard") + ggtitle("Mito percent") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol=3
)
ggsave(file=paste0(figures,"QC1.png"), p, width=3400, height=1400, unit="px")


## Plot cells with mitochondria metrics
p <- plotColData(unfiltered, x="sum", y="subsets_Mito_percent",
    colour_by="discard") + scale_x_log10()
ggsave(file=paste0(figures,"QC2.png"), p, width=1500, height=1500, unit="px")


# Normalization
sce.bone <- logNormCounts(sce.bone, size_factors = sce.bone$sum)
summary(sizeFactors(sce.bone))

# Feature selection
# Variance modeling
library(scran)
## Fit the trend to the variance
set.seed(1010010101)
dec.bone <- modelGeneVarByPoisson(sce.bone, 
    block=sce.bone$Donor, BPPARAM=bpp)
## Feature selection
top.bone <- getTopHVGs(dec.bone, n=5000)
saveRDS(top.bone, file=paste0(data, "top5000_hvgs.rds"))

## save variance plot
png(paste0(figures,"QC_variance_modeling.png"), width=1500, height=1500, unit="px", pointsize=18)
par(mfrow=c(4,2))
blocked.stats <- dec.bone$per.block
for (i in colnames(blocked.stats)) {
    current <- blocked.stats[[i]]
    plot(current$mean, current$total, main=i, pch=16, cex=0.5,
        xlab="Mean of log-expression", ylab="Variance of log-expression")
    curfit <- metadata(current)
    curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
}
dev.off()



# Data integration (this take a huge amount of time)
library(batchelor)
library(BiocNeighbors)
## run
set.seed(1010001)
merged.bone <- fastMNN(sce.bone, batch = sce.bone$Donor, subset.row = top.bone,
     BSPARAM=BiocSingular::RandomParam(deferred = TRUE), 
     BNPARAM=AnnoyParam(),
     BPPARAM=bpp)

reducedDim(sce.bone, 'MNN') <- reducedDim(merged.bone, 'corrected')

# We use the percentage of variance lost as a diagnostic measure
metadata(merged.bone)$merge.info$lost.var


# Dimensionality reduction (take a huge amount of time)
set.seed(01010100)
sce.bone <- runUMAP(sce.bone, dimred="MNN",
    external_neighbors=TRUE, 
    BNPARAM=AnnoyParam(),
    BPPARAM=bpp,
    n_threads=bpnworkers(bpp))

# Clustering
library(bluster)
## do the clustering
set.seed(1000)
colLabels(sce.bone) <- clusterRows(reducedDim(sce.bone, "MNN"),
    TwoStepParam(KmeansParam(centers=1000), NNGraphParam(k=5)))
## Number of cells in each cluster
table(colLabels(sce.bone))
## Number of cells in each cluster per Donor
tab <- table(Cluster=colLabels(sce.bone), Donor=sce.bone$Donor)

## Plot heatmap
library(pheatmap)
p <- pheatmap(log10(tab+10), color=viridis::viridis(100))
ggsave(file=paste0(figures,"pheatmap.png"), p, width=1500, height=1500, unit="px")


# TODO: add scrambling option in scater's plotting functions.
scrambled <- sample(ncol(sce.bone))

p <- gridExtra::grid.arrange(
    plotUMAP(sce.bone, colour_by="label", text_by="label"),
    plotUMAP(sce.bone[,scrambled], colour_by="Donor")
)
ggsave(file=paste0(figures,"UMAP.png"), p, width=1500, height=2500, unit="px")

# Save processed data
saveRDS(sce.bone, file=paste0(data, "sce_bone_processed.rds"))

# Differential Expression
markers.bone <- findMarkers(sce.bone, block = sce.bone$Donor, 
    direction = 'up', lfc = 1, BPPARAM=bpp)


top.markers <- markers.bone[["2"]]
best <- top.markers[top.markers$Top <= 10,]
lfcs <- getMarkerEffects(best)
saveRDS(lfcs, file=paste0(data, "lfcs.rds"))
# Plot heatmap
library(pheatmap)
p <- pheatmap(lfcs, breaks=seq(-5, 5, length.out=101))
ggsave(file=paste0(figures,"top_markers_cluster_2_heatmap.png"), p, width=1500, height=1500, unit="px")


# Cell type classification
## Aggregate cluster expression
se.aggregated <- sumCountsAcrossCells(sce.bone, id=colLabels(sce.bone))

## Import reference human cell atlas data 
library(celldex)
hpc <- HumanPrimaryCellAtlasData()

## Running annotation
library(SingleR)
anno.single <- SingleR(se.aggregated, ref = hpc, labels = hpc$label.main,
    assay.type.test="sum")
anno.single

## Assign annotated labels to clusters
anno.single$labels
colData(sce.bone)$label
levels(colData(sce.bone)$label)

factor_data <- factor(
                    colData(sce.bone)$label, 
                    levels = levels(colData(sce.bone)$label), 
                    labels = anno.single$labels)
## Add labels to sce
colData(sce.bone)$celltype <- factor_data
## Plot celltype
p <- plotUMAP(sce.bone, colour_by="celltype", text_by="celltype")
ggsave(file=paste0(figures,"UMAP_celltype.png"), p, width=1500, height=1500, unit="px")
# Save processed data
saveRDS(sce.bone, file=paste0(data, "sce_bone_processed_annotated.rds"))

# Session info
# Get the session information
si <- sessionInfo()
# Extract attached packages
attached_pkgs <- si$otherPkgs
# Function to extract package names and versions
extract_pkg_info <- function(pkg_list) {
  data.frame(
    Tools = names(pkg_list),
    Version = sapply(pkg_list, function(x) x$Version),
    stringsAsFactors = FALSE
  )
}
# Convert attached packages to dataframes
attached_df <- extract_pkg_info(attached_pkgs)
# Print the dataframe
rownames(attached_df) <- NULL
attached_df 
## save session to file
saveRDS(attached_df, file=paste0(data, "sessionInfo.rds"))