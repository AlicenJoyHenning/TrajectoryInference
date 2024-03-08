# This R Script covers a workshop given by Bioconductor
# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html

# [1] Install packages 
BiocManager::install(c("slingshot", "SingleCellExperiment", "RColorBrewer",
                       "scales", "viridis", "UpSetR", "pheatmap", "msigdbr",
                       "fgsea", "knitr", "ggplot2", "gridExtra", "tradeSeq"))


library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(scales)
library(viridis)
library(UpSetR)
library(pheatmap)
library(msigdbr)
library(fgsea)
library(knitr)
library(ggplot2)
library(gridExtra)
library(tradeSeq)

# We are accessing the data through an R package 
install.packages("remotes")
remotes::install_github("kstreet13/bioc2020trajectories")
library(bioc2020trajectories)

# [2] Upstream workflow 

sce <- bioc2020trajectories::importRawData()
importRawData()
