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
get_marker_genes <- function(scExp, y, n.cores){

    if (missing(n.cores)) {
        # number of cores to be used, the default is to use all but one cores
        n.cores = detectCores() - 1
    }
   
    registerDoParallel(n.cores)

    ncells = y$N.cells#number of cells

    N.cluster = y$N.pred_cluster

    ##remove those unnecessary (no changes) genes
    cat("Preprocessing...\n")
    my = scExp
    
    my = my[apply(my,1,function(x) sd(x)!=0),]
    my <- t(scale(t(my)))

    mat = as.matrix(my)

    D = nrow(mat)#number of genes

    uy = y$unique_pred_clusters#unique
    if(max(as.numeric(uy)) > length(uy)){#in case we use a larger value than the number of unique elements
        y$pred_clusters  = match(y$pred_clusters, uy)#we use the index instead
    }

    dt <- data.frame(yy = t(mat), ig = y$pred_clusters)#gene expression; cluster

    cat("Looking for marker genes...\n")

    #    D = 100
    total <- D
    # create progress bar
    pb <- txtProgressBar(min = 0, max = total, style = 3)


    ######parallel
    ginfo = foreach(i = 1:D, .combine = c)%dopar%{
        r = rank(mat[i, ])#gene rank
   
        s = aggregate(r~dt$ig, FUN = mean)#average over one gene across cells
        icluster = which.max(s$r)#the cluster with the maximum value
        x = r[as.numeric(dt$ig) == icluster]
        y = r[as.numeric(dt$ig) != icluster]
        p = wilcox.test(x, y)$p.value
    
        pred <- prediction(r, as.numeric(dt$ig) == icluster)
        auc <- unlist(performance(pred, "auc")@y.values)
    
    
        tt = c(auc, icluster, p)
    
        # update progress bar
        setTxtProgressBar(pb, i)
    
        # cat(c(auc, icluster, p), "\n")
        return(tt)
    }

    ginfo = matrix(ginfo, nrow = 3)
    ###############


    close(pb)

    ginfo <- data.frame(matrix(unlist(ginfo), ncol = 3, byrow = T))
    colnames(ginfo) = c("auc", "icluster", "pvalue")
    rownames(ginfo) = rownames(mat)

    ginfo$pvalue = p.adjust(ginfo$pvalue, method = "holm")

    ###in case some clusters do not have any qualified marker genes, we need to lower the AUC threshold (0.85)
    dt = data.table(ginfo)
    dt0 = dt[, list(maxauc = max(auc), minpvalue = min(pvalue)), by = icluster]

    ##adjusted AUC and p-value
    adpvalue = max(0.01, max(dt0$minpvalue))
    adauc = min(0.8, min(dt0$maxauc))
    
    cat("The adjusted p-avalue is", adpvalue, "and the adjusted AUROC is", adauc, "\n")

    x1 = which(ginfo$pvalue < adpvalue)

    x2 = which(ginfo$auc > adauc)

    # x1 = which(ginfo$pvalue < 0.01)
    # 
    # x2 = which(ginfo$auc > 0.85)

    xx = intersect(x1, x2)

    sginfo = ginfo[xx, ]#selected genes with p-values and their clusters

    return(sginfo)
}

