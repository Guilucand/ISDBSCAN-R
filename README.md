# ISDBSCAN-R


## Summary 
Densitiy based clustering algorithm (ISDBSCAN) for single-cell RNA-seq

## Roadmap
This is an initial implementation, under active development.


## Installation

```R 3.5
#With R 3.5
install.packages("devtools")
library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("beachmat", version = "3.8")
BiocManager::install("SingleCellExperiment", version = "3.8")

install_github("InfOmics/ISDBSCAN-R")
```

```R 3.6
#With R 3.6
install.packages("devtools")
library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("beachmat", version = "3.9")
BiocManager::install("SingleCellExperiment", version = "3.9")

install_github("InfOmics/ISDBSCAN-R")
```

## Test data

```R 
library(ISDBSCAN)
library(ggplot2)

can383 = as.matrix(read.table(system.file("extdata", "can383.txt", package = "ISDBSCAN")))
res_list = ISDBSCAN(can383, k = 11, stratif = TRUE)

plotData = data.frame(cbind(can383, res_list$clusters))
plotData <- setNames(plotData, c("x","y","Clusters"))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
ggplot(plotData, aes(x = x, y = y), colors(distinct = TRUE)) + geom_point(aes(colour = factor(Clusters), group = Clusters))
```
Beside the dataset in example there are two other datasets:
    - can474.txt (use k = 11 for testing)
    - t4-normalized (use k = 12 for testing)
