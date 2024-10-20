# libraries ---------------------------------------------------------------
#Load libraries
library(tidyverse)
library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(BayesSpace)
library(patchwork)

# read in data ------------------------------------------------------------
# read in the processed data of the brain slices
list_brain <- readRDS(file = "out/object/list_brain_all_BayesSpace1000_q10.rds")

# generate an ad hoc plot
test <- list_brain$V01

Idents(test) <- "BayesSpace"
SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE, ncol = 4,alpha = 0.5) +
  plot_annotation(title = "BayesSpace")

# plot the full UMAP
DimPlot(test,label = T)

# filter only the one shown in the shortlist panel
test_subset <- subset(test,subset = BayesSpace %in% c(1,4,9,7))

DimPlot(test_subset)
ggsave("out/image/test_subset.pdf",width = 5,height = 5)

test@meta.data
