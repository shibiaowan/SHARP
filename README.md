# SHARP
**_SHARP_**: <b><u><i>S</i></u></b>ingle-cell RNA-seq <b><u><i>H</i></u></b>yper-fast and <b><u><i>A</i></u></b>ccurate processing via ensemble <b><u><i>R</i></u></b>andom <b><u><i>P</i></u></b>rojection.

# Introduction: 

SHARP is a bioinformatics tool to process and analyze single-cell RNA-seq (scRNA-seq) data. Algorithmically, SHARP is based on ensemble random projection (RP) and multi-layer meta-clustering which can well preserve cell-to-cell distance in reduced-dimensional space. Compared to other existing tools, it has the following advantages: 

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
It is recommended, although not mandatory, to specify what type of the input expression matrix. The data type of input expression matrix for SHARP can be fragments/reads per kilo base per million mapped reads (FPKM/RPKM), counts per million mapped reads (CPM), transcripts per million (TPM) or unique molecule identifiers (UMI), or other read counts. For consistency, FPKM/RPKM values are converted into TPM values and UMI values are converted into CPM values. Here is one example for UMI-based input expression matrix: 
```{r}
res = SHARP(scExp, exp.type = "UMI") # SHARP will converts UMI input data into CPM-based values
```
By default, SHARP treats the input as either TPM or CPM format and will not make the data-type conversion.


# Pre-processing:
Pre-processing for SHARP includes two parts: <b>removing all-zero genes</b> and <b>log-transform</b>. For the former part, for relatively small-size (e.g., 10,000 single cells) datasets, SHARP, by default, will remove those genes whose expressions are zero across all cells or whose expressions are missing (i.e., "NA") for some cells. SHARP also gives the option of whether removing or not to users:

```{r}
res = SHARP(scExp, prep = TRUE) # remove all-zero or expression-missing genes
```

For the log-transform (e.g., log(scExp + 1)), SHARP will do it for most datasets. However, SHARP will not do the log-transform in the following case: when the <a href="https://en.wikipedia.org/wiki/Silhouette_(clustering)" target="_blank">Silhouttee index</a> of the no-log-transform case is no less than a threshold (by default, 0.75) and larger than that of the log-transform case for a randomly subsampled subset (e.g., by default, 100 single cells) of the original data. This can guarantee the optimized clustering performance of various kinds of input expression data. We also give the option of doing the log-transform test to users:

```{r}
res = SHARP(scExp, logflag = FALSE) # do not check whether log-transform is done or not (for saving time); simply do the log-transform
```


# Number of Single Cells:

Depending the number of single cells in the input data, SHARP will automatically determine whether data-partitioning and similarity-based meta-clustering (sMetaC) are adopted or not (Note that if data partitioning is not adopted, sMetaC is neither necessary). For the former, a sub-function SHARP_small will be used, whereas for the latter, SHARP_large will be used. Experiments suggest that these two strategies may affect the clustering performance, especially for small-size (<1,000 single cells) datasets, but may have little effect on relatively large-size datasets (>5,000 single cells). However, data partitioning is undoubtedly able to accelerate the clustering process for SHARP. By default, when the number of single cells is fewer than 5,000, SHARP will not adopt the data-partitioning and sMetaC strategies. You can change the threshold as follows:

```{r}
res = SHARP(scExp, base.ncells = 2000) # when the number of single cells is larger than 2000, data partitioning and sMetaC are adopted
```

By default, the number of cells in each partition is 2000. You can also change the partition cell number:

```{r}
res = SHARP(scExp, base.ncells = 2000, partition.ncells = 1000) # when data partitioning and sMetaC are adopted, the number of cells for each partition is set to 1,000
```


# Number of Clusters:

By default, SHARP will automatically determine the optimal number of clusters by integrating Silhouttee index, Calinski-Harabasz index and hierarchical heights. However, you can also pre-define the final number of clusters, the minimum number of clusters checked by SHARP, or the maximum number of clusters checked by SHARP, respectively:

```{r}
res = SHARP(scExp, N.cluster = 8) # predefine the final cluster number is 8

res = SHARP(scExp, minN.cluster = 8) # let SHARP find the optimal cluster number from the minimum 8 and a maximum number determined by SHARP

res = SHARP(scExp, maxN.cluster = 8) # let SHARP find the optimal cluster number from the minimum number determined by SHARP and the maximum 8

```

SHARP adopts hierarchical clustering as the basic clustering method in individual RP clustering and weighted-base meta-clustering (wMetaC). You can predefine the number of clusters independently for each of the two cases:

```{r}
res = SHARP(scExp, indN.cluster = 8) # predefine the cluster number for individual RP clustering (whose cell size may be the partition cell size) as 8

res = SHARP(scExp, enpN.cluster = 8) # predefine the cluster number for wMetaC as 8
```

# Number of Reduced Dimension:

SHARP uses RP as a dimension-reduction method, which can well preserve the cell-to-cell distance in the reduced dimension. To what extent the dimension can be reduced is an interesting question. By default, SHARP determines the number of reduced dimension as a function of the number of cells (e.g., ceiling(log2(ncells)/(0.2^2))). You can also determine it by your own:

```{r}
res = SHARP(scExp, reduced.ndim = 500) # predefine the reduced dimension as 500
```

# Number of Random-Projection Applications:

To make the performance robust, SHARP by default adopts 15 runs of RPs (or ensemble size) for small-size datasets, whereas 5 runs of RPs for large-size datasets. You may also want to change it for either speeding the clustering or yielding more robust performance as follows:

```{r}
res = SHARP(scExp, ensize.K = 7) # predefine the reduced dimension as 500
```


# Deterministic Results:

SHARP is based ensemble random projection which will produce robust yet stochastic clustering results. To achieve reproducible and deterministic results as follows:

```{r}
res = SHARP(scExp, rN.seed = 10) # achieve reproducible results
```

# Visualization:

SHARP can also be used for single-cell data visualization, which runs much faster and probably better than the tSNE package. This is because SHARP uses the ensemble of the ensemble dimension-reduced feature matrix and the clustering-result induced matrix, which have significantly lower dimension than the original one and at the same time incorporates the global information from the clustering results. SHARP visualization can be achieved as follows:

```{r}
res = SHARP(scExp) # clustering
visualization_SHARP(res, label) # label is the reference (or predefined) clustering label
```

# Marker Genes:

SHARP is capable of identifying marker genes.

```{r}
res = SHARP(scExp) # clustering
sginfo = get_marker_genes(scExp, res) # detect marker genes
```

# Multi-Core Processing:

By default, SHARP is configured in parallel computing using multiple cores, i.e., using (n-1) cores, where n is the number of cores of the host computer. You can choose the number of cores to be used:

```{r}
res = SHARP(scExp, n.cores = 1) # running in a single core
```


# Others:

# Processing 1.3 Million Single Cells:



Briefly speaking, SHARP integrates the following algorithms and techniques to achieve both efficient and effective performance: (1) divide-and-conquer strategy; (2) random-projection based dimension reduction; (3) ensemble learning for both meta-clustering of results from different random projections and of those from different partitions of large-scale datasets; (4) weighted ensemble clustering for tackling difficult-to-cluster single cells and (5) similarity-based meta-clustering for combining results of different mutually-exclusive single-cell groups.

# Citation:

Shibiao Wan, Junil Kim and Kyoung Jae Won. SHARP: Single-Cell RNA-Seq Hyper-Fast and Accurate Processing via Ensemble Random Projection, submitted, 2018.
