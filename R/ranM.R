#' Prepare a random matrix
#'
#' This function is to prepare a random matrix to be used in parallel computation, so that the same random matrix will be applied to different partitions/groups of datasets in the same application of random projection.
#' 
#' @param scdata input single-cell expression matrix
#' @param p the dimension to be reduced to
#'
#' @import Matrix
#'
#' @export
ranM <- function(scdata, p){#for random projection; note scdata: m*n, m is the feature dimensions, n is the sample number; p is the reduced dimension
  m = nrow(scdata)#the number of features
  n = ncol(scdata)#number of samples/cells
#   set.seed(123)#a flag to fix the random number
  s = sqrt(m)#according to the paper "Very Sparse Random Projection"
  x0 = sample(c(sqrt(s),0, -sqrt(s)), size= m*p, replace=TRUE, prob=c(1/(2*s), 1 - 1/s, 1/(2*s)))
#   return(x)
  x <- Matrix(x0, nrow = m, byrow = TRUE, sparse = TRUE)#reshape a vector to a sparse matrix
 
  return(x)#the same format, feature*sample
}