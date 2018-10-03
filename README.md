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

#`install.packages("devtools")`#if you have not installed the package "devtools"

```{r}
library(devtools)
install_github("shibiaowan/SHARP")
```


# Quick Start: 

Suppose your input scRNA-seq expression matrix is "scExp"
`library(SHARP)`

`res = SHARP(scExp)`


# Example: 
An example data (i.e., scExp_tpm.RData) has been loaded together with SHARP package. Get to know some information of the example data:

`dim(scExp_tpm)`#the number of genes, number of cells

`scExp_tpm[1:5,1:5]`#typical values of the expression matrix

SHARP provides choices on multiple configurations per users' requests:

1. By default, SHARP automatically determines all of the parameters and/or configurations
`results_small = SHARP(scExp_tpm)`

2. You can predefine the number of clusters
`results_small = SHARP(scExp_tpm, N.cluster = 8)`#we predefine the number of clusters as 8

3. If you 
try SHARP for small-size datasets (by default, "small-size dataset" means a dataset containing less than 2000 single cells)

`results_small = SHARP(scExp_tpm)`

#try SHARP for large-size datasets (suppose you think that a dataset containing 300+ single cells is already very large for your local computational resources)

`results_large = SHARP(scExp_tpm, base.ncells = 300)`


Briefly speaking, SHARP integrates the following algorithms and techniques to achieve both efficient and effective performance: (1) divide-and-conquer strategy; (2) random-projection based dimension reduction; (3) ensemble learning for both meta-clustering of results from different random projections and of those from different partitions of large-scale datasets; (4) weighted ensemble clustering for tackling difficult-to-cluster single cells and (5) similarity-based meta-clustering for combining results of different mutually-exclusive single-cell groups.

# Citation:

Shibiao Wan, Junil Kim and Kyoung Jae Won. SHARP: Single-Cell RNA-Seq Hyper-Fast and Accurate Processing via Ensemble Random Projection, submitted, 2018.
