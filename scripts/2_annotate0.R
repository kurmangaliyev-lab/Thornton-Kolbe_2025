library(tidyverse)
library(Seurat)
library(cowplot)

rm(list = ls())
setwd("~/projects1/projJC1/analysis1/")

# 1. Load & Prep
data0 <- readRDS("data/data2.rds")
markers0 <- readRDS("data/markers2.rds")

# 2. Visual System Atlas (V1.1)
dataOL0 <- readRDS("/work1/yerbol/atlas0/atlas_V1.1/data_V1.1a.rds")
dataOL0 <- subset(dataOL0, class1 == "N" & time == "48h")

## reorder clusters
Idents(dataOL0) <- dataOL0$type2
order0 <- intersect(levels(dataOL0), paste0("N", 0:200))
order0 <- c(setdiff(levels(dataOL0), order0), order0)
dataOL0$type2 <- factor(dataOL0$type2, levels = order0)
Idents(dataOL0) <- dataOL0$type2

# 3. Average Expression
## fru+
mat0 <- as.matrix(AverageExpression(data0)[["RNA"]])
expr0 <- mat0 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "cluster2", values_to = "expr") %>%
  as.data.frame()

## atlas
matOL0 <- as.matrix(AverageExpression(dataOL0)[["RNA"]])
exprOL0 <- matOL0 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "type2", values_to = "expr") %>%
  as.data.frame()

# 4. Match fru+ to visual neurons
# 4.1 Features
features1 <- markers0 %>%
  filter(pct.1 > 0.75, avg_log2FC >= 2, p_val_adj < 0.01) %>%
  group_by(cluster2) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1:10) %>%
  pull(gene) %>%
  unique() %>%
  intersect(rownames(matOL0))

# 4.2 Correlation
cor1 <- cor(log1p(mat0[features1,]), log1p(matOL0[features1,]), method = "pearson") %>%
  as.data.frame() %>%
  rownames_to_column(var = "cluster2") %>%
  pivot_longer(cols = -cluster2, names_to = "type2", values_to = "r") %>%
  mutate(cluster2 = factor(cluster2, levels(data0)), type2 = factor(type2, levels(dataOL0)))

cor2 <- cor1 %>%
  group_by(cluster2) %>%
  arrange(desc(r)) %>%
  mutate(rank1 = row_number()) %>%
  ungroup() %>%
  group_by(type2) %>%
  arrange(desc(r)) %>%
  mutate(rank2 = row_number()) %>%
  ungroup() %>%
  as.data.frame()

# 4.3 Rename
cor2 %>% filter(rank1 == 1, rank2 == 1)
rename1 <- cor2 %>% filter(r >= 0.9) %>% droplevels() %>% arrange(type2)
rename2 <- setNames(as.character(rename1$type2), rename1$cluster2)
rename2 <- c("X3"="KCgm","X35"="KCgd", "X123"="KCab", rename2)

data0@meta.data <- data0@meta.data %>% 
  mutate(type2 = recode_factor(cluster2, !!!rename2)) %>%
  mutate(type2 = factor(type2, c(rename2, setdiff(levels(data0), rename2))))
Idents(data0) <- data0$type2
expr0 <- expr0 %>% 
  mutate(type2 = recode_factor(cluster2, !!!rename2)) %>%
  mutate(type2 = factor(type2, levels(data0$type2))) %>%
  relocate(expr, .after = type2)

# 5. Plots
DimPlot(data0, label = T, group.by = "type2") + NoLegend()

# 5.1 Correlation
theme1 <- theme(
  axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5), 
  axis.text.y = element_text(size = 6), axis.title = element_blank(), legend.position = "bottom"
)
cor1 %>%
  filter(cluster2 %in% names(rename2)) %>%
  ggplot(aes(x = type2, y = cluster2, fill = r)) +
  geom_tile(color = "grey20") +
  scale_fill_distiller(palette = "Spectral") +
  theme1
ggsave(filename = "results/cor1.pdf", width = 12, height = 4)

cor1 %>%
  filter(type2 %in% rename2) %>%
  ggplot(aes(y = type2, x = cluster2, fill = r)) +
  geom_tile(color = "grey20") +
  scale_fill_distiller(palette = "Spectral") +
  theme1
ggsave(filename = "results/cor2.pdf", width = 12, height = 4)

# 5.2 Marker genes
features2 <- markers0 %>%
  mutate(type2 = recode_factor(cluster2, !!!rename2)) %>%
  filter(type2 %in% setdiff(rename2, c("KCgm", "KCgd", "KCab")), pct.1 > 0.75, avg_log2FC >= 2, p_val_adj < 0.01) %>%
  group_by(type2) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(1:5) %>%
  pull(gene) %>%
  unique() %>%
  intersect(rownames(dataOL0))
features2 <- unique(c("fru","ey","dac","Mef2","prt","sNPF","ab", "norpA","GstD1","DIP-kappa", features2))

heatmap1 <- function(pdata0){
  ggplot(pdata0, aes(x = gene, y = type2, fill = log1p(expr2))) +
    geom_tile(color = "grey20") +
    scale_fill_distiller(palette = "RdYlBu", direction = -1) +
    scale_x_discrete(limits = features2) +
    theme(
      axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 6),
      axis.title = element_blank(), legend.position = "bottom", legend.title = element_blank()
    )
}

expr0 %>%
  filter(type2 %in% rename2, gene %in% features2) %>%
  mutate(type2 = factor(type2, rename2), expr2 = ifelse(expr > 20, 20, expr)) %>%
  heatmap1()
ggsave(filename = "results/markers1.pdf", width = 12, height = 4)

exprOL0 %>%
  filter(type2 %in% rename2, gene %in% features2) %>%
  mutate(type2 = factor(type2, rename2), expr2 = ifelse(expr > 20, 20, expr)) %>%
  heatmap1()
ggsave(filename = "results/markers2.pdf", width = 12, height = 4)

# 7. Save
saveRDS(data0, file = "data/data3.rds")
saveRDS(expr0, file = "data/expr3.rds")
