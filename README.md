# SHARP

<b><u><i><font size="6"><span style="color:red">S</span></font></i></u></b>ingle-cell RNA-seq <b><u><i>H</i></u></b>yper-fast and <b><u><i>A</i></u></b>ccurate processing via ensemble <b><u><i>R</i></u></b>andom <b><u><i>P</i></u></b>rojection

### News:

* <b>Mar-01-2021</b>: The version 1.1.0 of SHARP has been released. The new version SHARP addressed the installation problem arising from the "clues" package removed from CRAN.
* <b>Apr-22-2020</b>: The SHARP_unlimited function has been updated for better processing ultra-large-scale single-cell data.
* <b>Feb-20-2020</b>: An example of detailing the installation time and configurations has been added.
* <b>Jun-28-2019</b>: The documentations of SHARP have been significantly enriched to include details of parameter settings for SHARP.
* <b>Apr-04-2019</b>: The flashClust function has been added for fast implementation of hierarchical clustering.

### Table of Contents

[Introduction](https://github.com/shibiaowan/SHARP/blob/master/README.md#introduction)

[Installation](https://github.com/shibiaowan/SHARP/blob/master/README.md#installation)

[Quick Start](https://github.com/shibiaowan/SHARP/blob/master/README.md#quick-start)

[More Details](https://github.com/shibiaowan/SHARP#more-details)<details><summary>...</summary>
* [Expression Type](https://github.com/shibiaowan/SHARP#expression-type)
* [Pre-Processing](https://github.com/shibiaowan/SHARP#pre-processing)
* [Number of Single Cells](https://github.com/shibiaowan/SHARP#number-of-single-cells)
* [Number of Clusters](https://github.com/shibiaowan/SHARP#number-of-clusters)

* [Number of Reduced Dimension](https://github.com/shibiaowan/SHARP#number-of-reduced-dimension)

* [Number of Random-Projection Applications](https://github.com/shibiaowan/SHARP#number-of-random-projection-applications)

* [Reproducible Results](https://github.com/shibiaowan/SHARP#reproducible-results)

* [Multi-Core Processing](https://github.com/shibiaowan/SHARP#multi-core-processing)

* [Performance Evaluation](https://github.com/shibiaowan/SHARP#performance-evaluation)

* [Running Multiple Times of SHARP](https://github.com/shibiaowan/SHARP#running-multiple-times-of-sharp)

* [hclust vs flashClust](https://github.com/shibiaowan/SHARP#hclust-vs-flashclust)

* [Others](https://github.com/shibiaowan/SHARP#others)
</details>

[Visualization](https://github.com/shibiaowan/SHARP#visualization)

[Marker Genes](https://github.com/shibiaowan/SHARP#marker-genes)

[Processing 1.3 Million Single Cells](https://github.com/shibiaowan/SHARP/blob/master/README.md#processing-13-million-single-cells)

[Citation](https://github.com/shibiaowan/SHARP/blob/master/README.md#citation)

[FAQs](https://github.com/shibiaowan/SHARP#faqs)

[Bug Report](https://github.com/shibiaowan/SHARP/blob/master/README.md#bug-report)


# Introduction: 

SHARP, short for <b><u><i>S</i></u></b>ingle-cell RNA-seq <b><u><i>H</i></u></b>yper-fast and <b><u><i>A</i></u></b>ccurate processing via ensemble <b><u><i>R</i></u></b>andom <b><u><i>P</i></u></b>rojection, is a bioinformatics tool to process and analyze single-cell RNA-seq (scRNA-seq) data. Algorithmically, SHARP is based on ensemble random projection (RP) and multi-layer meta-clustering which can well preserve cell-to-cell distance in reduced-dimensional space. Compared to other existing tools, it has the following advantages: 

1. <b>scalable</b> to processing <b>10 million single cells</b>; 
2. <b>hyper-fast</b>; 
3. <b>accurate</b>;
4. <b>robust</b>; 
5. <b>R-based</b>, which can be readily used and combined with other R-based scRNA-seq data analysis tools; and 
6. <b>multi-functional</b>, including clustering, dimension reduction, fast single-cell data visualization, marker gene identification, etc.


# Installation:

```{r}
install.packages("devtools")#only if you have not installed the package "devtools"

library(devtools)
install_github("shibiaowan/SHARP")

```
The installation process takes around 47 seconds for a typical computer like with configuration of Intel Core i5-7300U @ 2.60GHz with 64 bit OS and a hard disk of 500 GB

# Quick Start: 

Load the library:
```{r}
library(SHARP)
```
An example data (i.e., scExp_tpm.RData) has also been loaded and we can check some basic information about it as follows:
```{r}
scExp = scExp_tpm #a TPM-based single-cell pancreas data from Wang et al.
dim(scExp) #check the numbers of genes and cells

scExp[1:5,1:5] #check typical values of the expression matrix
```
Run SHARP:
```{r}
res = SHARP(scExp)

#the last several lines should be something as follows:
#
##Analysis complete!
##-----------------------------------------------------------------------
##Total running time: 0.16285088 minutes
```
The running time for the example is less than 11 seconds for a computer with configuration of Intel Core i5-7300U @ 2.60GHz with 64 bit OS and a hard disk of 500 GB. Then, you can check the perfornance using the ARI (adjusted Rand index) by comparing your prediction results with the given reference clustering label of the example:
```{r}
ARI(label, res) 
# Rand        HA        MA        FM   Jaccard
#0.9675143 0.9201827 0.9206001 0.9435416 0.8918298
```
The variable "label" has already been loaded automatically when starting using SHARP.

If you want to yield reproducible results (e.g., the results above) by SHARP, you can set a fixed number for the parameter rN.seed:
```{r}
rN = 2103; 
res = SHARP(scExp_tpm, rN.seed=rN); 
ARI(label, res)
```

# More Details:

## Expression Type:
It is recommended, although not mandatory, to specify what type of the input expression matrix. The data type of input expression matrix for SHARP can be fragments/reads per kilo base per million mapped reads (FPKM/RPKM), counts per million mapped reads (CPM), transcripts per million (TPM) or unique molecule identifiers (UMI), or other read counts. For consistency, FPKM/RPKM values are converted into TPM values and UMI values are converted into CPM values. Here is one example for UMI-based input expression matrix: 
```{r}
res = SHARP(scExp, exp.type = "UMI") # SHARP will converts UMI input data into CPM-based values
```
By default, SHARP treats the input as either TPM or CPM format and will not make the data-type conversion.


## Pre-Processing:
Pre-processing for SHARP includes two parts: <b>removing all-zero genes</b> and <b>log-transform</b>. For the former part, for relatively small-size (e.g., 10,000 single cells) datasets, SHARP, by default, will remove those genes whose expressions are zero across all cells or whose expressions are missing (i.e., "NA") for some cells. SHARP also gives the option of whether removing or not to users:

```{r}
res = SHARP(scExp, prep = TRUE) # remove all-zero or expression-missing genes
```

For the log-transform (e.g., log(scExp + 1)), SHARP will do it for most datasets. However, SHARP will not do the log-transform in the following case: when the <a href="https://en.wikipedia.org/wiki/Silhouette_(clustering)" target="_blank">Silhouttee index</a> of the no-log-transform case is no less than a threshold (by default, 0.75) and larger than that of the log-transform case for a randomly subsampled subset (e.g., by default, 100 single cells) of the original data. This can guarantee the optimized clustering performance of various kinds of input expression data. We also give the option of doing the log-transform test to users:

```{r}
res = SHARP(scExp, logflag = FALSE) # do not check whether log-transform is done or not (for saving time); simply do the log-transform
```


## Number of Single Cells:

Depending on the number of single cells in the input data, SHARP will automatically determine whether data-partitioning and similarity-based meta-clustering (sMetaC) are adopted or not (Note that if data partitioning is not adopted, sMetaC is neither necessary). For the former, a sub-function SHARP_small will be used, whereas for the latter, SHARP_large will be used. Experiments suggest that these two strategies may somewhat affect the clustering performance, especially for small-size (<1,000 single cells) datasets, but may have little effect on relatively large-size datasets (>5,000 single cells). However, data partitioning is undoubtedly able to significantly accelerate the clustering process for SHARP, especially for large-scale datasets. By default, when the number of single cells is fewer than 5,000, SHARP will not adopt the data-partitioning and sMetaC strategies. You can change the threshold by "base.ncells" as follows:

```{r}
res = SHARP(scExp, base.ncells = 2000) # when the number of single cells is larger than 2000, data partitioning and sMetaC are adopted
```

By default, the number of cells in each partition is 2000. You can also change the partition cell number by "partition.ncells":

```{r}
res = SHARP(scExp, base.ncells = 2000, partition.ncells = 1000) # when data partitioning and sMetaC are adopted, the number of cells for each partition is set to 1,000
```


## Number of Clusters:

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

## Number of Reduced Dimension:

SHARP uses RP as a dimension-reduction method, which can well preserve the cell-to-cell distance in the reduced dimension. To what extent the dimension can be reduced is an interesting question. By default, SHARP determines the number of reduced dimension as a function of the number of cells (e.g., ceiling(log2(ncells)/(0.2^2))). You can also determine it by your own:

```{r}
res = SHARP(scExp, reduced.ndim = 500) # predefine the reduced dimension as 500
```

## Number of Random-Projection Applications:

To make the performance robust, SHARP by default adopts 15 runs of RPs (or ensemble size) for small-size datasets, whereas 5 runs of RPs for large-size datasets. You may also want to change it for either speeding the clustering or yielding more robust performance as follows:

```{r}
res = SHARP(scExp, ensize.K = 7) # predefine the reduced dimension as 500
```


## Reproducible Results:

SHARP is based ensemble random projection which will produce robust yet stochastic clustering results. To achieve reproducible and deterministic results as follows:

```{r}
res = SHARP(scExp, rN.seed = 10) # achieve reproducible results
```


## Multi-Core Processing:

By default, SHARP is configured in parallel computing using multiple cores, i.e., using (n-1) cores, where n is the number of cores of the host computer. You can choose the number of cores to be used as follows:

```{r}
res = SHARP(scExp) # running with multi-cores, where (n-1) cores are used and n is the total number of CPU cores


res = SHARP(scExp, n.cores = 1) # running in a single core
```


## Performance Evaluation:

SHARP evalulates the clustering performance based on Adjusted Rand Index (ARI) and running time. The running time will be given as soon as SHARP finishes clustering. For the ARI, SHARP calculates as follows:

```{r}
res = SHARP(scExp)
ARI(q, res) # q is the reference (or ground-truth) clustering label
```
SHARP will give 5 different ARI-based metrics, in which "HA" (Hubert and Arabie's ARI) is corresponding to the ARI we often use. 


## Running Multiple Times of SHARP:

SHARP will produce robust yet stochastic clustering results. To evaluate the clustering performance, running multiple times of SHARP is sometimes necessary. SHARP provides this kind of functions as follows:
```{r}
allres = run_Mtimes_SHARP(scExp, Mtimes = 10) # run 10 times of SHARP
```

## hclust vs flashClust:

The <a href="https://cran.r-project.org/web/packages/flashClust/index.html" target="_blank">flashClust</a> package is a fast implementation of hierarchical clustering. It is claimed to run faster than the traditional hierarchical clustering package (i.e., hclust) for large-scale datasets. Therefore, for users' convenience, we have also embedded flashClust into SHARP. We added the option (i.e., flashmark) in the SHARP package to allow users to select whether flashClust (i.e., flashmark = TRUE) or the traditional hclust (i.e., flashmark = FALSE) is used as the base clustering method for SHARP. SHARP uses flashClust as the base clustering method as follows: 

```{r}
res = SHARP(scExp, flashmark = TRUE) # use the flashClust as the base clustering method
```
By default, flashmark = FALSE, i.e., hclust is used as the base clustering method.

## Others:

SHARP uses the hierarchical clustering (hclust) as the basic clustering method. The parameters related to hclust are also customizable, including hierarchical method (e.g., "ward.D", "complete", "single", etc).

```{r}
res = SHARP(scExp, hmethod = "ward.D") # using "ward.D" as the hierarchical clustering method
```

SHARP integrates Silhouttee index, CH index and hierarchical heights to determine the optimal number of clusters. The threshold for using CH index instead of Silhouttee index and the threshold of the hierarchical height difference for cutting trees are also available to adjust.

```{r}
res = SHARP(scExp, sil.thre = 0.5) # when the avearge Silhouttee index is less than 0.5 (by default, 0.35), CH index will be used to optimize the cluster number

res = SHARP(scExp, height.Ntimes = 3) # when the descending-ordered hierarchical height is 3 times larger than the next immediate height, we cut the tree 
```

If the number of single cells is very large, even the dimension-reduced feature matrix can be very large. In case that scatter-plot visualization is not necessary, we do not need to store the dimension-reduced feature matrix in the clustering results. The following option can achieve this goal:

```{r}
res = SHARP(scExp, forview = FALSE) # do not save the feature matrix
```


# Visualization:

SHARP can also be used for single-cell data visualization, which runs much faster and probably better than the tSNE package. This is because SHARP uses the ensemble of the ensemble dimension-reduced feature matrix and the clustering-result induced matrix, which have significantly lower dimension than the original one and at the same time incorporates the global information from the clustering results. SHARP visualization can be achieved as follows:

```{r}
res = SHARP(scExp) # clustering
visualization_SHARP(res, label) # label is the reference (or predefined) clustering label
```

# Marker Genes:

SHARP is capable of identifying marker genes. It uses a differential approach which checks the statistical significance of cluster-specific genes vs genes in other clusters. Then, it uses cluster-specific sparsity ratio to exclude those less significant genes. A heatmap of the several top cluster-specific marker genes can be drawn.

```{r}
res = SHARP(scExp) # clustering
sginfo = get_marker_genes(scExp, res) # detect marker genes

sortmarker = plot_markers(sginfo) # plot marker gene heatmap
```


# Processing 1.3 Million Single Cells:

One of the unique contributions for SHARP is that it can process a dataset with 1.3 million single cells. To the best of our knowledge, SHARP is the first R package which can process and analyze more than one million single cells. Existing scRNA-seq tools can't tackle it because the 1.3 million scRNA-seq data matrix can not be directly loaded into R due to the lack of 64-bit integers support in R. On the contrary, SHARP can handle it. Because SHARP adopted the divide-and-conquer strategy to analyze huge datasets, we first divided the 1,306,127 scRNA-seq data into 26 smaller-size blocks of data, i.e., each block with 50,000 single cells except the last one with 56,127 single cells. Next, SHARP analyzed each block of data, whose results were integrated by our proposed algorithm sMetaC. If we save the 1.3 million scRNA-seq data in a list format, we can use the following command to deal with:

```{r}
res = SHARP_unlimited(scExp) # dealing with 1.3 million single cells which are saved as a list of 26 matrices
```

For memory-efficient processing, we can also save each block of data into an RDS file (RDS files are suggested because they are compact and their sizes are smaller) and each time we only process one file. In this case, we just need to provide the directory of those files

```{r}
ndinfo = list()#three elements: the directory of those files to save the 1.3 million single cells, number of cells and number of genes
ndinfo$dir = "tmp/million_cells/"
ndinfo$ncells = 1306127
ndinfo$ngenes = 27998
res = SHARP_unlimited3(ndinfo) # dealing with 1.3 million single cells which are saved as a list of 26 matrices
```

After clustering, we can identify the marker genes as follows:

```{r}
sginfo = get_marker_genes_unlimited2(gdinfo, res) # detect marker genes
```


# Citation:

Shibiao Wan, Junil Kim and Kyoung Jae Won. SHARP: Single-Cell RNA-Seq Hyper-Fast and Accurate Processing via Ensemble Random Projection, <b><i>Genome Research</i></b>, 2020, vol. 30, pp. 205-213. ([link](https://genome.cshlp.org/content/30/2/205)) ([Preprint](https://www.biorxiv.org/content/early/2018/11/04/461640))

# FAQs:

__Q__: When installing SHARP, I came across a problem stating that 'Error in {: task 1 failed-"there is no "get_CH""'. What should I do?

__A__: Please install the latest version of SHARP. The latest version (i.e., v1.1.0) of SHARP has addressed the problem from the package "clues" which was removed from CRAN, so it should be cleared off the problem aforementioned.

<br/>

__Q__: When installing SHARP, I came across a problem stating that "ERROR: dependency 'clues' is not available for package 'SHARP'". What should I do?

__A__: Please install the latest version of SHARP. The latest version (i.e., v1.1.0) of SHARP has addressed the problem from the package "clues" which was removed from CRAN, so it should be cleared off the problem aforementioned.

<br/>

__Q__: When installing SHARP, I came across a problem stating that "clang: error: unsupported option '-fopenmp'" and/or "ERROR: compilation failed for package ‘data.table’". What should I do?

__A__: Please select "no" when prompted about a question "Do you want to install from sources the package which needs compilation? (Yes/no/cancel)".


# Bug Report:

If you find any bugs or problems, or you have any comments on SHARP, please don't hesitate to contact the software package maintaner Shibiao Wan at shibiaowan.work@gmail.com.


# Copyright

Copyright © 2018-2020 [Shibiao Wan](https://sites.google.com/site/shibiaowan/)




# License 

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

GNU GENERAL PUBLIC LICENSE  
Version 3, 29 June 2007

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
