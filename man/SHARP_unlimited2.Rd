% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SHARP_unlimited2.R
\name{SHARP_unlimited2}
\alias{SHARP_unlimited2}
\title{Run SHARP for clustering single-cell RNA data whose size (e.g., dim > 2,147,483,647) is beyond R to handle.}
\usage{
SHARP_unlimited2(
  scExp,
  ensize.K,
  reduced.ndim,
  partition.ncells,
  hmethod,
  N.cluster = NULL,
  enpN.cluster = NULL,
  indN.cluster = NULL,
  minN.cluster,
  maxN.cluster,
  sil.thre,
  height.Ntimes,
  logflag,
  n.cores,
  forview = TRUE,
  rN.seed
)
}
\arguments{
\item{scExp}{a list of input single-cell expression matrices, where each element represents a partition of the huge-size single-cell RNA-seq data matrix.}

\item{viewflag}{a logic to indicate whether to save the ensemble random-projection feature vectors for further visualization. The default is TRUE.}

\item{...}{other parameters similar to those used in SHARP() function. Please refer to SHARP() for details.}
}
\value{
a list containing the SHARP clustering results, distribution of the clustering results, the predicted optimal number of clusters, time SHARP consumes for clustering, some intermediate results including clustering results by each random-projection based hierarchical clustering and other related statstical information including number of cells, genes, reduced dimensions and number of applications of random projection.
}
\description{
Because the current R lacks 64-bit integers support, it is recommended to divided the huge-size data into several smaller partitions first.
}
\details{
This is a complementary SHARP function to deal with the case where R itself can not handle a matrix whose non-zero element is over 2,147,483,647. The huge-size matrix is first divided into several partitions which are susequently saved into a list. SHARP_unlimited() can directly deal with this list.
}
\examples{
enresults = SHARP_unlimited(scExp)

}
\author{
Shibiao Wan <shibiao@pennmedicine.upenn.edu>
}
