library(tidyverse)
library(Seurat)
library(UCell)
library(cowplot)

rm(list = ls())
setwd("~/projects1/projJC1/analysis1/")

# 1. Load and Prep: Palmateer et al. 2023
data0 <- readRDS("data/GSE160370_FullDataSet_Merge.rds")
data0 <- UpdateSeuratObject(data0)

meta0 <- data0@meta.data
counts0 <- GetAssayData(data0, assay = "RNA", layer = "counts")
rename0 <- c(
  "CG42313"="side-II", "CG34113"="side-III", "CG14372"="side-IV", "CG34371"="side-V", "CG34114"="side-VI", "CG12950"="side-VII", 
  "CG12484"="side-VIII", "CG31814"="DIP-kappa", "CG45781"="DIP-lambda", "CG10824"="cDIP", "CG15630"="fipi", "CG1149"="MstProx"
)
rownames(counts0) <- recode(rownames(counts0), !!!rename0)

data0 <- CreateSeuratObject(counts = counts0, project = "Palmateer_2023")
data0 <- AddMetaData(data0, metadata = meta0 %>% select(perc_mt = percent.mt, replicate, sex, type1 = FullAnnotation))

Idents(data0) <- data0$replicate
saveRDS(data0, "data/data1.rds")

# 2. QC/Filter
data0 <- readRDS("data/data1.rds")
VlnPlot(data0, c("nCount_RNA")) + scale_y_continuous(trans = "log10", breaks = 500*2^c(0:5)) + geom_hline(yintercept = 1000)
data0 <- subset(data0, nCount_RNA > 1000)
dim(data0)

# 3. Clustering
dims0 <- 1:100

# 3.1 Cluster 1
data0 <- data0 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c("orig.ident", "nCount_RNA")) %>%
  RunPCA(npcs = 200) %>%
  RunTSNE(dims = dims0) %>%
  #RunUMAP(dims = dims0) %>%
  FindNeighbors(dims = dims0) %>%
  FindClusters(resolution = 2, cluster.name = "cluster1")
data0$cluster1 <- factor(paste0("X", data0$cluster1), levels = paste0("X", levels(data0)))

# 3.2 Cluster 2
data0 <- data0 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  AddModuleScore_UCell(
    features = list(
      oxp = c("ATPsynbeta", "ATPsyngamma"), 
      hsp = c("Hsp67Ba", "Hsp67Bc"),
      arg = c("Hr38", "sr", "CG14186")
    )
  ) %>%
  ScaleData(vars.to.regress = c("orig.ident", "nCount_RNA", "oxp_UCell", "hsp_UCell", "arg_UCell")) %>%
  RunPCA(npcs = 200) %>%
  RunTSNE(dims = dims0) %>%
  #RunUMAP(dims = dims0) %>%
  FindNeighbors(dims = dims0) %>%
  FindClusters(resolution = 2, cluster.name = "cluster2")
data0$cluster2 <- factor(paste0("X", data0$cluster2), levels = paste0("X", levels(data0)))

## Plot Clusters
plot_grid(
  DimPlot(data0, group.by = "cluster1", shuffle = T, label = T, label.size = 3) + NoLegend() + NoAxes(),
  DimPlot(data0, group.by = "cluster2", shuffle = T, label = T, label.size = 3) + NoLegend() + NoAxes(),
  plot_grid(
    FeaturePlot(data0, features = "Hr38", label = F, min.cutoff = "q10", max.cutoff = "q90") + NoAxes(),
    FeaturePlot(data0, features = "ATPsynbeta", label = F, min.cutoff = "q10", max.cutoff = "q90") + NoAxes(),
    FeaturePlot(data0, features = "Hsp67Bc", label = F, min.cutoff = "q10", max.cutoff = "q90") + NoAxes(),
    FeaturePlot(data0, features = "alphaTub84B", label = F, min.cutoff = "q10", max.cutoff = "q90") + NoAxes()
  ),
  nrow = 1
)

# 4. Save
Idents(data0) <- data0$cluster2
saveRDS(data0, "data/data2.rds")

# 5. Markers
markers0 <- FindAllMarkers(data0, min.pct = 0.25, logfc.threshold = 1)
markers0 <- markers0 %>% select(cluster2 = cluster, gene, everything())
saveRDS(markers0, "data/markers2.rds")
