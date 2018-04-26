#' Run multiple-times SHARP for single-cell RNA data clustering
#'
#' This function is to run multiple times of SHARP for evaluating SHARP in clustering of single-cell RNA-Seq data
#'
#' @param scExp input single-cell expression matrix
#' @param ensize.K number of applications of random projection for ensemble
#' @param reduced.ndim the dimension to be reduced to
#' @param base.ncells a base threshold of number of cells. When the number of cells of a dataset is smaller than this threshold, we use SHARP_small function; otherwise, we use SHARP_large.
#' @param partition.ncells number of cells for each partition when using SHARP_large
#' @param n.cores number of cores to be used. The default is (n-1) cores, where n is the number of cores in your local computer or server.
#' @param Mtimes number of times to run SHARP
#'
#' @examples
#' enresults = Run_Mtimes_SHARP(scExp, ensize.K, reduced.ndim, partition.ncells, base.ncells, n.cores, Mtimes)
#'
#' @import foreach
#'
#' @import parallel
#'
#' @import doParallel
#'
#' @export
Run_Mtimes_SHARP <- function(scExp, ensize.K, reduced.ndim, partition.ncells, base.ncells, 
    n.cores, Mtimes) {
    rm(list = ls())  #remove all
    
    ########## Input matrix files###########
    
    # Files = c('Wang/Pancreas_Julia.published.dataTpms_selected.txt') Files =
    # c('Enge/tpms4MATLAB_selected.txt')#you need to make transpose of it because it
    # is feature * samples Files = c('Goolam/Goolam_expression_log10.txt') Files =
    # c('Park/Park_CPM.rds')
    A = Mtimes
    K = ensize.K
    
    # eps = c(0.3, 0.25, 0.2, 0.15, 0.1, 0.05)
    allresults <- vector("list", length = length(K))  #declare a matrix of lists
    
    for (k in K) {
        # allresults = list() for(i in 1:length(eps)){
        info = vector("list", length = A)
        kname = paste("enSize_", k, sep = "")
        # dim = ceiling(log2(2282)/(eps[i]^2)) kname = paste('dim_', dim, sep = '')
        for (j in 1:A) {
            # times of using SHARP
            print("=========================================================================", 
                quote = FALSE)
            print(paste("SHARP: Ensemble Size-", k, ", Run-", j, sep = ""), quote = FALSE)
            cat("=========================================================================\n")
            cat("SHARP: Ensemble Size-", k, ", Run-", j, "\n")
            # enresults = getallresults_large(Files, K = k)#using the default Dim
            enresults = SHARP(scExp, k, reduced.ndim, partition.ncells, base.ncells, 
                n.cores)
            jname = paste("Run_", j, sep = "")
            info[[jname]] = enresults
        }
        allresults[[kname]] = info
        # } # dset = enresults$dataname#just use the last one to retrieve the dataset
        # name
    }
    
    
    # print('Done!', quote = FALSE)
    cat("Done!\n")
    return(allresults)
}
