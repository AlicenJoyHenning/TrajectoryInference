# This R Script covers a workshop given by Bioconductor
# https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html

# [1] Install packages ####
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
library(patchwork)

# We are accessing the data through an R package 
install.packages("remotes")
remotes::install_github("kstreet13/bioc2020trajectories")
library(bioc2020trajectories)

# [2] Upstream workflow ####
getwd()
setwd("C:/Users/alice/OneDrive/Documents/GitHub/bioc2020trajectories")
# sce <- bioc2020trajectories::importRawData()

# [3] Workflow Process ####
# Load data
data("sce", package = "bioc2020trajectories")

# After integration, the first question is of differential topology: 
# Q: are there very different topologies btw the 2 conditions 

# creates a new variable 'shuffle' of integer value
shuffle <- sample(ncol(sce)) 
layout(matrix(1:2, nrow = 1))
par(mar = c(4.5, 4, 1, 1))

plot(reducedDims(sce)$UMAP[shuffle, ],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = alpha(c(1:2)[factor(colData(sce)$pheno$treatment_id)][shuffle], alpha = .5))
legend("topright", pch = 16, col = 1:2, bty = "n",
       legend = levels(factor(colData(sce)$pheno$treatment_id)))
plot(reducedDims(sce)$UMAP[shuffle, ], 
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(3, 4)[factor(colData(sce)$pheno$spatial_id)][shuffle], alpha = .5))
legend("topright", pch = 16, col = 3:4, bty = "n", 
       legend = levels(factor(colData(sce)$pheno$spatial_id)))
layout(1)
par(mar = c(5, 4, 4, 2) + .1)

# Next, looking at the balance of the data 
help("imbalance_score", "bioc2020trajectories")
# Compute an imbalance score to show whether nearby cells have the same condition of not

scores <- bioc2020trajectories::imbalance_score(
  rd = reducedDims(sce)$UMAP, # reduced dimension matrix of the cells
  cl = colData(sce)$pheno$treatment_id, # vector of conditions 
  k = 20, # number of neighbours to consider (default 10)
  smooth = 40 # smoothing parameter where lower values mean more smoothing
  )

grad <- viridis::plasma(10, begin = 0, end = 1)
names(grad) <- levels(cut(scores$scaled_scores, breaks = 10))
plot(reducedDims(sce)$UMAP, col = grad[cut(scores$scaled_scores, breaks = 10)],
     asp = 1, pch = 16, xlab = "UMAP_1", ylab = "UMAP_2", cex = .8)
legend("topleft", legend = names(grad), col = grad, pch = 16, bty = "n", cex = 2 / 3)
# we see that start & end position are the same regardless of treatment condition 
# meaning we will only be fitting one trajectory 


# Trajectory Analysis 
# Use slingshot for trajectory inference with the cells position (inner/outer)
# acting as the cluster identifier. 

sce <- slingshot(sce, reducedDim = "UMAP", 
                 clusterLabels = colData(sce)$pheno$spatial_id, 
                 start.clus = 'inner', approx_points = 150)

grad <- viridis::viridis(10, begin = 0, end = 1)
plot(reducedDims(sce)$UMAP, col = grad[cut(scores$scaled_scores, breaks = 10)], pch = 16, asp = 1)
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')




# Differential progression 

ks.test(slingPseudotime(sce)[colData(sce)$pheno$treatment_id == "Mock", 1],
        slingPseudotime(sce)[colData(sce)$pheno$treatment_id == "TGFB", 1])

# Differential expression 

set.seed(3)
icMat <- evaluateK(counts = as.matrix(assays(sce)$counts),
                   pseudotime = colData(sce)$slingshot$pseudotime,
                   cellWeights = colData(sce)$slingshot$cellWeights.V1,
                   conditions = factor(colData(sce)$pheno$treatment_id),
                   nGenes = 300,
                   k = 3:7)

# ERROR : Either provide the slingshot object using the sds argument, or provide pseudotime and cell-level weights manually using pseudotime and cellWeights arguments.

set.seed(3)
sce <- fitGAM(sce,
              conditions = factor(colData(sce)$pheno$treatment_id),
              nknots = 5)
# | 1 % ~02h 44m 43s   takes long time 
mean(rowData(sce)$tradeSeq$converged)

# Assess DE along pseudotime
rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))

assocRes <- rowData(sce)$assocRes
mockGenes <- rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionMock, "fdr") <= 0.05)
]
tgfbGenes <- rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionTGFB, "fdr") <= 0.05)
]

length(mockGenes)# tells us the number of DEGs in the mock dataset
length(tgfbGenes)# tells us the number of DEGs in the treatment condition
UpSetR::upset(fromList(list(mock = mockGenes, tgfb = tgfbGenes)))

# Visualization of DEGs 

yhatSmooth <- predictSmooth(sce, gene = mockGenes, nPoints = 50, tidy = FALSE)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE, 
                       show_row_names = FALSE, 
                       show_colnames = FALSE)

# Cluster genes according to their average expression pattern 
cl <- sort(cutree(heatSmooth$tree_row, k = 6))
table(cl)
## cl 
## 1    2    3    4    5    6
## 

conditions <- colData(sce)$pheno$treatment_id 
pt1 <- colData(sce)$slingshot$pseudotime







