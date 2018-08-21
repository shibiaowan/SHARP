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
#' enresults = SHARP_unlimited(scExp)
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
SHARP_unlimited <- function(scExp, viewflag = TRUE, ...) {
    # timing
    start_time <- Sys.time()  #we exclude the time for loading the input matrix
    
    title = "scRNA-Seq Clustering"
    
    if (missing(scExp)) {#scRNA-seq expression matrix
        stop("No expression data is provided!")
    }
    
    if(class(scExp) != "list"){
        if(class(scExp) == "matrix"){
            warning("SHARP is used instead of SHARP_unlimited because the input is a matrix!")
            enresults = SHARP(scExp)
            return(enresults)
        }
        stop("The input should be a LIST of partitioned scRNA-seq expression matrices!")
    }else if(class(scExp) == "list" && length(scExp) == 1){
        warning("SHARP is used instead of SHARP_unlimited because the length of the input is 1!")
        enresults = SHARP(scExp[[1]])
        return(enresults)
    }
    
#     if (missing(n.cores)) {
#         # number of cores to be used, the default is to use all but one cores
#         n.cores = detectCores() - 1
#     }
#     # registerDoMC(n.cores)
#     registerDoParallel(n.cores)

    nnp = length(scExp)#number of partitions/blocks
    nng = vector("numeric", nnp)
    nnc = vector("numeric", nnp)
    
#     ncells = sum(sapply(1:length(scExp), function(i) scExp[[i]]@Dim[2]))#number of cells
    ncells = sum(sapply(1:length(scExp), function(i) dim(scExp[[i]])[2]))#number of cells
    p = ceiling(log2(ncells)/(0.2^2))#the reduced dimension
    
    fColor = character(0)
    E1 = numeric(0)
    folds = numeric(0)
    y = list()
    for(i in 1:nnp){
        cat("=======================================================================\nProcessing Partition", i, "of the scRNA-seq data...\n")
        nng[i] = nrow(scExp[[i]]) 
        nnc[i] = ncol(scExp[[i]])
        
        folds = c(folds, rep(i, nnc[i]))#folds
        
        mat = scExp[[i]]
        
#         mat = as.matrix(scExp[[i]])
        cat("=======================================================================\n")
        
        y[[i]] = SHARP(mat, reduced.ndim = p, prep = FALSE, logflag = FALSE, ...)
        
        fColor = c(fColor, paste(y[[i]]$pred_clusters, "s", i, sep = ""))
        E1 = rbind(E1, y[[i]]$viE)
#         if (!missing(exp.type)) {#expression type
#             if (exp.type != "CPM" && exp.type != "TPM" && exp.type != "UMI") {#if the type is neither CPM nor TPM
#                 scExp = t(t(scExp)/colSums(scExp)) * 1e+06#we do normalization
#             }
#         }
    }
    
    cat("Length of fColor:", length(fColor), "\n")
    cat("Number of unique elements of fColor:", length(unique(fColor)), "\n")
    cat("Dim of E1:", dim(E1), "\n")
#     print(y[[1]]$paras)
    
    cat("=======================================================================\n")
    cat("Meta-clustering for combining blockwise results...\n")
    stmp = sMetaC(fColor, E1, folds, y[[1]]$paras$hmethod, y[[1]]$paras$N.cluster, y[[1]]$paras$minN.cluster, 
            y[[1]]$paras$maxN.cluster, y[[1]]$paras$sil.thre, y[[1]]$paras$height.Ntimes)
    finalrowColor = stmp$finalColor
    
    
    
    
    
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
    
    enresults = list()
    enresults$pred_clusters = finalrowColor
    enresults$unique_pred_clusters = uf
    enresults$distr_pred_clusters = tn1
    enresults$N.pred_clusters = luf
    enresults$N.cells = sum(nnc)
    enresults$N.genes = nng[1]
    enresults$reduced.dim = y[[1]]$reduced.ndim
    enresults$ensize.K = y[[1]]$ensize.K
    if(viewflag){
        enresults$viE = E1
        
#         x0 = matrix(0, ncells, luf)
#         for(i in 1:ncells){
#             x0[i, finalrowColor[i] == uf] = 1
#         }
        x0 = sparseMatrix(i = 1:ncells, j = as.numeric(finalrowColor), x = 1, dims = c(ncells, luf))
        enresults$x0 = x0
    }
    enresults$time = t
    enresults$paras = y[[1]]$paras
    
    return(enresults)
}
