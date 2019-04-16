# ISDBSCAN-R


## Summary 
Densitiy based clustering algorithm (ISDBSCAN) for single-cell RNA-seq

## Roadmap
This is an initial implementation, under active development.


## Installation

```R
install.packages("devtools")
library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("beachmat", version = "3.8")
BiocManager::install("SingleCellExperiment", version = "3.8")

install_github("InfOmics/ISDBSCAN-R")
```
