# SHARP
**_SHARP_**: <b><u><i>S</i></u></b>ingle-cell RNA-Seq <b><u><i>H</i></u></b>yper-fast and <b><u><i>A</i></u></b>ccurate clustering via ensemble <b><u><i>R</i></u></b>andom <b><u><i>P</i></u></b>rojection.

# Installation: 

`**library(devtools)**`

`**install_github("shibiaowan/SHARP")**`

# Quick Start: 

`**library(SHARP)**`#an example data (i.e., scExp_tpm.RData) has been loaded together with SHARP package

#get to know the number of genes and number of cells of the example data

`**dim(scExp_tpm)**`

#try SHARP for small-size datasets

`**results_small = SHARP(scExp_tpm)**`

#try SHARP for large-size datasets

`**results_large = SHARP(scExp_tpm, base.ncells = 300)**`

# Introduction: 

SHARP is a bioinformatics tool to process and analyze single-cell RNA-Seq (scRNA) data  in a hyper-fast and accurate way. It can process a scRNA dataset of 47,000+ single cells within several minutes with accurate clustering performance. 

Briefly speaking, SHARP integrates the following algorithms and techniques to achieve both efficient and effective performance: (1) divide-and-conquer strategy; (2) random-projection based dimension reduction; (3) ensemble learning for both meta-clustering of results from different random projections and of those from different partitions of large-scale datasets; (4) weighted ensemble clustering for tackling difficult-to-cluster single cells and (5) similarity-based meta-clustering for combining results of different mutually-exclusive single-cell groups.

# Citation:

Shibiao Wan, Junil Kim and Kyoung Jae Won. SHARP: Single-Cell RNA-Seq Hyper-Fast and Accurate Clustering via Ensemble Random Projection, submitted, 2018.
