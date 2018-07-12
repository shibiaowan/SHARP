#' Obtain a random matrix and its related information
#'
#' This function is to generate a random matrix according to a sparse-version Achlioptas distribution (see our paper for details), given the input matrix and the dimension we want to reduce to.
#'
#' @param scdata input single-cell expression matrix
#' @param p the dimension to be reduced to
#'
#' @examples
#' newE = RPmat(E, p)
#'
#' @import Matrix
#'
#' @export
RPmat <- function(scdata, p, seedn) {
    # for random projection; note scdata: m*n, m is the feature dimensions, n is the
    # sample number; p is the reduced dimension
    m = nrow(scdata)  #the number of features
    n = ncol(scdata)  #number of samples/cells
    # set.seed(123)#a flag to fix the random number
    s = sqrt(m)  #according to the paper 'Very Sparse Random Projection'
    if (seedn%%1 == 0) {
        # integer
        set.seed(seedn)
        x0 = sample(c(sqrt(s), 0, -sqrt(s)), size = m * p, replace = TRUE, prob = c(1/(2 * 
            s), 1 - 1/s, 1/(2 * s)))
    } else {
        x0 = sample(c(sqrt(s), 0, -sqrt(s)), size = m * p, replace = TRUE, prob = c(1/(2 * 
            s), 1 - 1/s, 1/(2 * s)))
    }
    # return(x)
    x <- Matrix(x0, nrow = m, byrow = TRUE, sparse = TRUE)  #reshape a vector to a sparse matrix
    projmat = 1/sqrt(p) * t(x) %*% scdata  #the projected matrix by random projection;matrix multiplication %*%
    # rmat = ranmat(m, p)#the random matrix
    
    # if(n <10000){ x <- matrix(x0, nrow = m, byrow = TRUE)#reshape a vector to a
    # matrix # scmat = data.matrix(scdata)#convert the data frame to a matrix projmat
    # = 1/sqrt(p)*t(x)%*%scdata#the projected matrix by random projection;matrix
    # multiplication %*% }else{ x2 <- Matrix(x0, nrow = m, byrow = TRUE, sparse =
    # TRUE)#reshape a vector to a matrix projmat = 1/sqrt(p)*t(x2)%*%scdata#the
    # projected matrix by random projection;matrix multiplication %*% }
    
    # projmat = t(projmat)
    list1 = list()
    list1$R = x
    list1$projmat = projmat
    return(list1)  #the same format, feature*sample
}
