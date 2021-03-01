#' Run SHARP for clustering single-cell RNA data whose size (e.g., dim > 2,147,483,647) is beyond R to handle.
#'
#' Because the current R lacks 64-bit integers support, it is recommended to divided the huge-size data into several smaller partitions first.
#'
#' @param scExp a list of input single-cell expression matrices, where each element represents a partition of the huge-size single-cell RNA-seq data matrix.
#'
#' @param viewflag a logic to indicate whether to save the ensemble random-projection feature vectors for further visualization. The default is TRUE.
#'
#' @param ... other parameters similar to those used in SHARP() function. Please refer to SHARP() for details.
#'
#' @details This is a complementary SHARP function to deal with the case where R itself can not handle a matrix whose non-zero element is over 2,147,483,647. The huge-size matrix is first divided into several partitions which are susequently saved into a list. SHARP_unlimited() can directly deal with this list.
#'
#' @return a list containing the SHARP clustering results, distribution of the clustering results, the predicted optimal number of clusters, time SHARP consumes for clustering, some intermediate results including clustering results by each random-projection based hierarchical clustering and other related statstical information including number of cells, genes, reduced dimensions and number of applications of random projection.
#'
#' @examples
#' enresults = SHARP_unlimited3(ndinfo)
#'
#' @author Shibiao Wan <shibiao@pennmedicine.upenn.edu>
#'
#' @import foreach
#'
#' @import parallel
#'
#' @import doParallel
#'
#' @import Matrix
#'
#' @export
SHARP_unlimited3 <- function(ndinfo, viewflag = TRUE, n.cores, ensize.K, rN.seed, N.cluster = NULL, ...) {
    # timing
    start_time <- Sys.time()  #we exclude the time for loading the input matrix
    
    title = "scRNA-Seq Clustering"
    
     if (missing(ndinfo)) {#scRNA-seq expression matrix
        stop("No expression data is provided!")
    }
    
    scExp_dir = ndinfo$dir#folder containing the partitions of the single cells datasets
    ncells = ndinfo$ncells#number of cells given
    ngenes = ndinfo$ngenes#number of genes
    if(!dir.exists(scExp_dir)){#the folder to store the results
	  stop(paste0(scExp_dir, " should be a folder storing several partitions of single-cell datasets!\n"))
        } 
        s = scExp_dir
        if(substr(s, nchar(s), nchar(s)) == "/"){#the last character of the dir is "/"
            s = substr(s, 1, nchar(s)-1)#remove the last "/"
        }
        scExp_dir = s
 
    if (missing(n.cores)) {
        # number of cores to be used, the default is to use all but one cores
        n.cores = detectCores() - 1
    }
    # registerDoMC(n.cores)
    registerDoParallel(n.cores)

    
    allfiles = list.files(scExp_dir, full.names = TRUE)#include the full path
    x1 = as.numeric(gsub("\\D*([0-9]+).*$", "\\1", allfiles))#
    allfiles = allfiles[order(x1)]#reorganize the file order according to the file name with digits
    len = length(allfiles)
    nnp = len
    
    
#     nnp = length(scExp)#number of partitions/blocks
    nng = vector("numeric", nnp)
    nnc = vector("numeric", nnp)
    
#     ncells = sum(sapply(1:length(scExp), function(i) scExp[[i]]@Dim[2]))#number of cells
#     ncells = sum(sapply(1:length(scExp), function(i) dim(scExp[[i]])[2]))#number of cells
    p = ceiling(log2(ncells)/(0.2^2))#the reduced dimension
    
    if (!missing(rN.seed)) {#random seed for reproducible results
        # seed number for reproducible results
        if (!is.numeric(rN.seed)) {
            stop("The rN.seed should be a numeric!")
        }
        if (rN.seed%%1 != 0) {
            stop("The rN.seed should be an integer!")
        }
    } else {
        rN.seed = 0.5  #use a non-integer to represent random results
    }
    
    if (missing(ensize.K)) {
                ensize.K = 5
    }
            
    rM = foreach(k = 1:ensize.K) %dopar% {###to use consistent random matrices across different partitions
            if (rN.seed == 0.5) {
                ranM2(ngenes, p, 0.5)  #for random results
            } else {
                ranM2(ngenes, p, 50 + rN.seed + k)  #for reproducible result
            }
    }
#     stopImplicitCluster()
        
    ifColor = vector("list", length = nnp)
    iE1 = vector("list", length = nnp)
    ifolds = vector("list", length = nnp)
    y = list()
    for(i in 1:nnp){
        cat("=======================================================================\nProcessing Partition", i, "of the scRNA-seq data...\n")
        mat = readRDS(allfiles[i])
        cat("Process", allfiles[i], "\n")
        nng[i] = nrow(mat) 
        nnc[i] = ncol(mat)
    
        
#         mat = as.matrix(scExp[[i]])
        cat("=======================================================================\n")
        
        y[[i]] = SHARP(mat, reduced.ndim = p, prep = FALSE, n.cores = n.cores, rM = rM, ensize.K = ensize.K, rN.seed = rN.seed, ...)
        
#         folds = c(folds, rep(i, nnc[i]))#folds
#         fColor = c(fColor, paste(y[[i]]$pred_clusters, "s", i, sep = ""))
#         E1 = rbind(E1, y[[i]]$viE)
        
        ifolds[[i]] = rep(i, nnc[i])#folds
        ifColor[[i]] = paste(y[[i]]$pred_clusters, "s", i, sep = "")
        iE1[[i]] = t(y[[i]]$viE)###for ease of ordering; by default, R concatenate elements by columns; but here the feature dim (column) is the same for each element of the list
        
        rm(mat)
        gc()
#         if (!missing(exp.type)) {#expression type
#             if (exp.type != "CPM" && exp.type != "TPM" && exp.type != "UMI") {#if the type is neither CPM nor TPM
#                 scExp = t(t(scExp)/colSums(scExp)) * 1e+06#we do normalization
#             }
#         }
    }
    
#     folds = unlist(ifolds)#convert a list to a vector
#     fColor = unlist(ifColor)
#     E1 = t(matrix(unlist(iE1), ncol = ncells))
    
    folds = foreach(i = 1:nnp, .combine = c)%dopar%{
        ifolds[[i]]
    }
    
    fColor = foreach(i = 1:nnp, .combine = c)%dopar%{
        ifColor[[i]]
    }
    E1 = foreach(i = 1:nnp, .combine = cbind)%dopar%{
        iE1[[i]]
    }
    E1 = t(E1)
    
    cat("=======================================================================\n")
    cat("Length of fColor:", length(fColor), "\n")
    cat("Number of unique elements of fColor:", length(unique(fColor)), "\n")
    cat("Dim of E1:", dim(E1), "\n")
#     print(y[[1]]$paras)
    
    cat("=======================================================================\n")
    cat("Meta-clustering for combining blockwise results...\n")
    stmp = sMetaC(fColor, E1, folds, y[[1]]$paras$hmethod, N.cluster, y[[1]]$paras$minN.cluster, 
            y[[1]]$paras$maxN.cluster, y[[1]]$paras$sil.thre, y[[1]]$paras$height.Ntimes)
    finalrowColor = stmp$finalColor
    
    
    print(table(finalrowColor))
    if(!is.numeric(N.cluster) && ncells > 1e4){#remove those extremely small clusters
        cat("Adjust clusters with very small number of cells...\n")
        xt = table(finalrowColor)
        s = names(which(xt < 10))#combine those clusters whose numbers are less than 10
        if(length(s) != 0){
             xi = which(!is.na(match(finalrowColor, as.numeric(s))))
#         s1 = length(unique(finalrowColor)) - length(s) + 1
            finalrowColor[xi] = min(as.numeric(s))
        }
    }
    
    #reorganize the clustering in order according to the decreasing order of number of cells
    x = sort(table(finalrowColor), decreasing = TRUE) 
    map = setNames(as.character(1:length(x)), names(x))#the former replacing the latter 
    finalrowColor[] <- map[finalrowColor]
    finalrowColor = as.numeric(finalrowColor)
    
    
    end_time <- Sys.time()
    
    t <- difftime(end_time, start_time, units = "mins")  #difference time in minutes
    # print(t)
    
    cat("Analysis complete!\n")
    cat("=======================================================================\n")
    cat("Total running time:", t, "minutes\n")
    #################################### 
    uf = unique(finalrowColor)
    luf = length(uf)#length of unique pred-clusters
    
    tn = table(finalrowColor)
    tn1 = tn[order(as.numeric(names(tn)))]#order the distributions of the clustering results
    print(tn1)
    
    enresults = list()
    enresults$pred_clusters = finalrowColor
    enresults$unique_pred_clusters = sort(as.numeric(uf))
    enresults$distr_pred_clusters = tn1
    enresults$N.pred_clusters = luf
    enresults$N.cells = sum(nnc)
    enresults$N.genes = nng[1]
    enresults$reduced.dim = y[[1]]$reduced.ndim
    enresults$ensize.K = y[[1]]$ensize.K
    if(viewflag){
        if(ncells>1e5){#if the number cells is too large, we need to do further dim reduction
            kdim = 50
            if (rN.seed == 0.5) {
                z0 = ranM2(p, kdim, 0.5)  #for random results
            } else {
                z0 = ranM2(p, kdim, 50 + rN.seed + k)  #for reproducible result
            }
            
            enresults$viE = as.matrix(1/sqrt(kdim)*E1%*%z0)
        }else{#otherwise, we just used the original one
            enresults$viE = E1
        }

#         x0 = matrix(0, ncells, luf)
#         for(i in 1:ncells){
#             x0[i, finalrowColor[i] == uf] = 1
#         }
#         cat(length(1:ncells), "and", length(as.numeric(finalrowColor)), "\n")
        x0 = sparseMatrix(i = 1:ncells, j = as.numeric(finalrowColor), x = 1, dims = c(ncells, luf))
        enresults$x0 = x0
    }
    enresults$time = t
    enresults$paras = y[[1]]$paras
    
    
    stopImplicitCluster()
    
    return(enresults)
}
