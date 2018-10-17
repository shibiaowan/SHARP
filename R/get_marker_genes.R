#' Get marker genes for each cluster
#'
#' This is to identify marker genes for each cluster after SHARP clustering. It uses a method similar to that proposed in the SC3 package except two improvements: (1) it uses an adjusted threshold to select marker genes instead of using a hardthreshod (i.e., p-value < 0.01 and AUROC > 0.85), so that marker genes can be found for all clusters; (2) it uses a parallelization way to calculate the p-value and AUROC (areas under the curve of the Receiver Operating Characteristic) for each gene, thus much faster than SC3.
#'
#' @param scExp the original expression matrix
#'
#' @param y the clustering results after running SHARP.R
#'
#' @param n.cores   number of cores to be used. The default is (n-1) cores, where n is the number of cores in your local computer or server.
#'
#' @return a matrix where each row represents each selected marker gene and each column represnts one property of the gene, including its cluster, p-value and AUROC.
#'
#' @examples
#'
#' y = SHARP(scExp)
#' sginfo = get_marker_genes(scExp, y)
#'
#' @import ROCR
#'
#' @import doParallel
#'
#' @import data.table
#'
#' @export
get_marker_genes <- function(scExp, y, theta, auc, pvalue, ng, n.cores){

    # timing
    start_time <- Sys.time()  #we exclude the time for loading the input matrix
    if (missing(n.cores)) {
        # number of cores to be used, the default is to use all but one cores
        n.cores = detectCores() - 1
    }
   
    registerDoParallel(n.cores)

    if(missing(theta)){
        theta = 1e-4
    }
    
    if(missing(auc)){
        auc = 0.7
    }
    
    if(missing(pvalue)){
        pvalue = 0.01
    }
    
    if(missing(ng)){
        ng = 5
    }
   
    ncells = y$N.cells#number of cells

    N.cluster = y$N.pred_cluster

    ##remove those unnecessary (no changes) and duplicated genes
    cat("Preprocessing...\n")
    r = duplicated(rownames(scExp))#find any duplicated genes in the input expression
    dr = which(r)#duplicated indices
    if(length(dr)>0){
        warning(paste(length(dr), "duplicated genes are found and then are removed!"))
        scExp = scExp[!r, ]
    }
    
#     scExp = scExp[rowSums(scExp)!=0, ]#much faster than using standard deviation
#     scExp = scExp[apply(scExp,1,function(x) sd(x)!=0),]
    #remove those all-zero or all-same-expression genes
#     f0 = foreach(i = 1:nrow(scExp), .combine = c)%dopar%{
#         x = sd(scExp[i, ])
#         if(x!=0){
#             return(TRUE)
#         }else{
#             return(FALSE)
#         }
#     }
#     scExp = scExp[f0, ]
    
#     scExp <- t(scale(t(scExp)))#no need to normalize

#     scExp = as.matrix(scExp)

    D = nrow(scExp)#number of genes

    uy = y$unique_pred_clusters#unique
    if(max(as.numeric(uy)) > length(uy)){#in case we use a larger value than the number of unique elements
        y$pred_clusters  = match(y$pred_clusters, uy)#we use the index instead
    }

    label = as.numeric(y$pred_clusters)
#     dt <- data.frame(yy = t(scExp), ig = label)#gene expression; cluster
    dt <- data.frame(ig = label)#gene expression; cluster
    
    cat("Looking for marker genes...\n")

    #    D = 100
    total <- D
    # create progress bar
    pb <- txtProgressBar(min = 0, max = total, style = 3)

    gn0 = rownames(scExp)

    rr = min(ng, N.cluster)#number of top genes whose AUC to be checked
    ######parallel
#     ginfo = foreach(i = 1:D, .combine = c)%dopar%{
    g4 = foreach(i = 1:D, .combine = "rbind")%dopar%{
        dd = scExp[i, ]
        dp = length(which(dd!=0))/length(dd)
        
        if(dp > theta){
            r = rank(scExp[i, ])#gene rank
   
   
             s = aggregate(r~dt$ig, FUN = mean)#average over one gene across cells
#         cat("Get the cell index of the maximum expression...\n")
#             icluster = which.max(s$r)#the cluster with the maximum value
            
            iclusters = order(-s$r)[1:rr]#indices of the top rr average ranks
            auc_res0 <- foreach(icluster = iclusters, .combine = c)%do%{
                pred = prediction(r, as.numeric(dt$ig) == icluster)
                auc_res <- unlist(performance(pred, "auc")@y.values)
                return(auc_res)
            }
            
            ic = which.max(auc_res0)
            icluster = iclusters[ic]#optimal cluster ID
            auc0 = auc_res0[ic]#the maximum AUC
#             cat("AUC:", auc_res0, "--maxAUC:", auc0, "--icluster:", icluster, "\n")
            
            x1 = r[as.numeric(dt$ig) == icluster]
            x2 = r[as.numeric(dt$ig) != icluster]
            p = wilcox.test(x1, x2)$p.value
            
            
#             s = aggregate(r~dt$ig, FUN = mean)#average over one gene across cells
#             icluster = which.max(s$r)#the cluster with the maximum value
#             x = r[as.numeric(dt$ig) == icluster]
#             y = r[as.numeric(dt$ig) != icluster]
#             p = wilcox.test(x, y)$p.value
#     
#             pred <- prediction(r, as.numeric(dt$ig) == icluster)
#             auc0 <- unlist(performance(pred, "auc")@y.values)
        }else{
            auc0 = icluster = 0
            p = 1
        }
        
    
    
#         tt = c(auc, icluster, p)
        tt = c(auc0, icluster, p, dp)#character
    
        # update progress bar
        setTxtProgressBar(pb, i)
    
        # cat(c(auc, icluster, p), "\n")
        return(tt)
    }
    
#     gn1 = as.character(ginfo[,1])
#     ginfo = ginfo[,c(2:5)]
    
#     ginfo = matrix(ginfo, nrow = 4)
    ###############


    close(pb)

    rownames(g4) = gn0
    colnames(g4) = c("auc", "icluster", "pvalue", "sparsity")

    ginfo = data.frame(g4)
    
#     ginfo <- data.frame(matrix(unlist(ginfo), ncol = 4, byrow = T))
#     colnames(ginfo) = c("auc", "icluster", "pvalue", "sparsity")
#     rownames(ginfo) = rownames(scExp)

    nr = nrow(ginfo)
    cat("Number of genes checked:", nr, "\n")
    
    ginfo = ginfo[ginfo$sparsity > theta, ]#remove low-expressed genes
    
    ginfo = ginfo[!is.nan(ginfo$pvalue), ]#remove those corresponding to p-value == NaN
    
    ginfo$pvalue = p.adjust(ginfo$pvalue, method = "holm")

    
    gallinfo = ginfo
    colnames(gallinfo) = c("auc", "icluster", "pvalue", "sparsity")
    ###in case some clusters do not have any qualified marker genes, we need to lower the AUC threshold (0.85)
    dt = data.table(ginfo)
    dt0 = dt[, list(maxauc = max(auc), minpvalue = min(pvalue)), by = icluster]

    ##adjusted AUC and p-value
#     adpvalue = max(0.01, max(dt0$minpvalue))
#     adauc = min(0.85, min(dt0$maxauc))
    adpvalue = pvalue
    adauc = min(auc, min(dt0$maxauc))
    
    cat("The adjusted p-avalue is", adpvalue, "and the adjusted AUROC is", adauc, "\n")

    x1 = which(ginfo$pvalue < adpvalue)

    x2 = which(ginfo$auc > adauc)

    # x1 = which(ginfo$pvalue < 0.01)
    # 
    # x2 = which(ginfo$auc > 0.85)

    xx = intersect(x1, x2)

    sginfo = ginfo[xx, ]#selected genes with p-values and their clusters
    
    nsgene = nrow(sginfo)
    cat("Number of selected marker genes:", nsgene, "\n")
    #the corresponding expressions of marker genes
#     dd = foreach(i = 1:nsgene, .combine = "rbind")%dopar%{
#         unlist(sapply(1:len, function(x) scExp[[x]][rownames(sginfo[i, ]),]))
#     }
    
#     dd = scExp[rownames(sginfo), ]
#     dd = scExp[which(gn0 %in% rownames(sginfo)), ]
    
    dd = scExp[rownames(sginfo), ]#after removing the duplicated genes
    
    end_time <- Sys.time()
    
    t <- difftime(end_time, start_time, units = "mins")  #difference time in minutes
    cat("Running time for finding marker genes:", t ,"minutes\n")
    
    res = list()
    res$mginfo = sginfo
    res$mat = dd
    res$label = label
    res$gallinfo = gallinfo
    return(res)
#     return(sginfo)
}

