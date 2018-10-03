# SHARP
**_SHARP_**: <b><u><i>S</i></u></b>ingle-cell RNA-seq <b><u><i>H</i></u></b>yper-fast and <b><u><i>A</i></u></b>ccurate processing via ensemble <b><u><i>R</i></u></b>andom <b><u><i>P</i></u></b>rojection.

# Introduction: 

SHARP is a bioinformatics tool to process and analyze single-cell RNA-seq (scRNA-seq) data. Algorithmically, SHARP is based on ensemble random projection and multi-layer meta-clustering which can well preserve cell-to-cell distance in reduced-dimensional space. Compared to other existing tools, it has the following advantages: 

1. <b>scalable</b> to processing <b>1.3 million single cells</b>; 
2. <b>hyper-fast</b>; 
3. <b>accurate</b>;
4. <b>robust</b>; and 
5. <b>multi-functional</b>, including clustering, dimension reduction, fast single-cell data visualization, marker gene identification, etc.


# Installation:

```{r}
install.packages("devtools")#only if you have not installed the package "devtools"

library(devtools)
install_github("shibiaowan/SHARP")
```


# Quick Start: 

Load the library:
```{r}
library(SHARP)
```
An example data (i.e., scExp_tpm.RData) has also been loaded and we can check some basic information about it:
```{r}
scExp = scExp_tpm #a TPM-based single-cell pancreas data from Wang et al.
dim(scExp) #check the numbers of genes and cells

`scExp[1:5,1:5]` #check typical values of the expression matrix
```
Run SHARP:
```{r}
res = SHARP(scExp)
```


# Expression Type:

# Pre-processing:

# Number of Single Cells:

# Number of Clusters:

# Number of Reduced Dimension:

# Number of Random-Projection Applications:

# Number of Base Cells and Partition Cells:

# Visualization:

# Marker Genes:

# Multi-Core Processing:

# Others:

# Processing 1.3 Million Single Cells:

try SHARP for small-size datasets (by default, "small-size dataset" means a dataset containing less than 2000 single cells)

`results_small = SHARP(scExp_tpm)`

#try SHARP for large-size datasets (suppose you think that a dataset containing 300+ single cells is already very large for your local computational resources)

`results_large = SHARP(scExp_tpm, base.ncells = 300)`


Briefly speaking, SHARP integrates the following algorithms and techniques to achieve both efficient and effective performance: (1) divide-and-conquer strategy; (2) random-projection based dimension reduction; (3) ensemble learning for both meta-clustering of results from different random projections and of those from different partitions of large-scale datasets; (4) weighted ensemble clustering for tackling difficult-to-cluster single cells and (5) similarity-based meta-clustering for combining results of different mutually-exclusive single-cell groups.

# Citation:

Shibiao Wan, Junil Kim and Kyoung Jae Won. SHARP: Single-Cell RNA-Seq Hyper-Fast and Accurate Processing via Ensemble Random Projection, submitted, 2018.
