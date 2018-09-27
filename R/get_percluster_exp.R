#' Get cluster-specific expression matrix
#'
#' This is to get the cluster-specific expresion matrix for further detecting cluster-specific marker genes, given the original scRNA-seq data (in the list format) and the clustering results.
#'
#' @param scExp the list of original block-wise expression matrices
#'
#' @param y the clustering results by running SHARP
#'
#' @param n.cores   number of cores to be used. The default is (n-1) cores, where n is the number of cores in your local computer or server.
#'
#' @return a list of sparse matrices where each matrix corresponds to the expression matrix of each cluster
#'
#' @examples
#'
#' y = SHARP_unlimited(scExp)
#' smat = get_percluster_exp(scExp, y)
#'
#'
#' @import doParallel
#'
#' @import Matrix
#'
#' @export
get_percluster_exp <- function(scExp, y, n.cores){

    start_time <- Sys.time()  #we exclude the time for loading the input matrix
    
    l = length(scExp)
    
    if(missing(n.cores)){
        n.cores = min(l, detectCores()-1)
    }
    registerDoParallel(n.cores)
#     ind = which(label == ll)

#     ll = sort(unique(as.numeric(label)))
#     uc = length(ll)#number of unique clusters
    
    label = as.numeric(y$pred_clusters)
    uc = length(unique(label))
    lens = sapply(1:l, function(x){ncol(scExp[[x]])})#numbers of cells in each partition
    xx = sapply(1:l, function(x){sum(lens[1:x])})#threshold of numbers of cells; sum of previous parititions
    
    xx = c(0, xx)
    
    

    cat("Organize the matrix within each partition...\n")
    #return the i-th partition and the j-th column
    sE = foreach(k = 1: l)%dopar%{
        t1 = label[(xx[k]+1): xx[k + 1]]
        t2 = sort(unique(t1))
        se = sapply(1:length(t2), function(x){scExp[[k]][, which(t1 == t2[x]), drop = FALSE]})
        f = vector("list", length = uc)
        f[t2] = se#save cluster-specific matrix in the corresponding list matrix
        return(f)
    }
    
    cat("Get cluster-specific matrices...\n")
    smat = foreach(j = 1:uc)%do%{
        cat("Cluster", j, "\n")
        sk = foreach(k = 1:l, .combine = 'cfun', .multicombine=TRUE)%dopar%{
            sE[[k]][[j]]
        }
        return(sk)
    }
    
        
    stopImplicitCluster()
    
    end_time <- Sys.time()
    
    t <- difftime(end_time, start_time, units = "mins")  #difference time in minutes
    cat("Running time:", t ,"minutes\n")
    
    return(smat)
    
}

#' Combine results which may contain NULL or sparse matrices
#'
#' Because by default, combining of sparse matrices and NULL will yield NULL which is not consistent with our expectation, we need to remove NULL to do the combining the sparse matrices. This function acts as a self-crafted combination function for the foreach output.
#'
#' @param ... variable number of inputs which are either sparse matrices or NULL
cfun <- function(...){
        ina = list(...)
        len = length(ina)
        g = unlist(sapply(1:len, function(i){if(!is.null(ina[[i]])){return(i)}else{return(NULL)}}))
        if(length(g) != 0){
            newinda = ina[g]
            x = foreach(j = 1:length(newinda),.combine = cbind)%do%{newinda[[j]]}
        }else{
            x = NULL
        }
        return(x)
}



