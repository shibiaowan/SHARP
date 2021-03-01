#' Run multiple-times SHARP for single-cell RNA data clustering
#'
#' This function is to run multiple times of SHARP for evaluating SHARP in clustering of single-cell RNA-Seq data
#'
#' @param scExp input single-cell expression matrix
#' @param Mtimes number of times to run SHARP. By default, 10 times will be tried
#' @param Kset a set of ensemble sizes (ensize.K) to be tried. By default, only one element (i.e., ensize.K = 15) will be tried
#' @param ... other parameters of SHARP that have not been listed here
#'
#' @examples
#' enresults = run_Mtimes_SHARP(scExp, Mtimes = 10, Kset = 5, partition.ncells = 2000, n.cores = 1)
#'
#' @import foreach
#'
#' @import parallel
#'
#' @import doParallel
#'
#' @export
run_Mtimes_SHARP <- function(scExp, Mtimes = 10, Kset = 15, ...) {
    
    ########## Input matrix files###########
    
    # Files = c('Wang/Pancreas_Julia.published.dataTpms_selected.txt') Files =
    # c('Enge/tpms4MATLAB_selected.txt')#you need to make transpose of it because it
    # is feature * samples Files = c('Goolam/Goolam_expression_log10.txt') Files =
    # c('Park/Park_CPM.rds')
    A = Mtimes
    K = Kset
    
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
            enresults = SHARP(scExp, ensize.K = k, forview = FALSE, ...)
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
