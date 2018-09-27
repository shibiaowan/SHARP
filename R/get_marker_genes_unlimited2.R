#' Get marker genes of each cluster for huge-size single-cell datasets
#'
#' This is an version of finding marker genes for huge-size single-cell datasets (e.g., the size of the scRNA-seq data is so huge that the expression matrix can not be read into R; instead, a list of block-wise expression matrices is read into R.) This also corresponds to SHARP_unlimited(). Detailed specifications can be found in get_marker_genes().
#'
#' @param scExp the list of original block-wise expression matrices
#'
#' @param y the clustering results after running SHARP.R
#'
#' @param n.cores   number of cores to be used. The default is (n-1) cores, where n is the number of cores in your local computer or server.
#'
#' @return a matrix where each row represents each selected marker gene and each column represnts one property of the gene, including its cluster, p-value and AUROC.
#'
#' @examples
#'
#' y = SHARP_unlimited(scExp)
#' sginfo = get_marker_genes_unlimited2(scExp, y)
#'
#' @import ROCR
#'
#' @import doParallel
#'
#' @import data.table
#'
#' @export
get_marker_genes_unlimited2 <- function(gdinfo, y, theta, auc, pvalue, n.cores){

    start_time <- Sys.time()  #we exclude the time for loading the input matrix
    
    if (missing(n.cores)) {
        # number of cores to be used, the default is to use all but one cores
        n.cores = detectCores() - 1
    }
   
    registerDoParallel(n.cores)
    
    if(missing(theta)){
        theta = 1e-5
    }
    
    if(missing(auc)){
        auc = 0.85
    }
    
    if(missing(pvalue)){
        pvalue = 0.05
    }

    ncells = y$N.cells#number of cells

    N.cluster = y$N.pred_cluster

#     ##remove those unnecessary (no changes) genes
    cat("Preprocessing...\n")
#     len = length(scExp)
    D = y$N.genes#number of genes
    cat("Number of all genes:", D, "\n")

    s = gdinfo
    if(substr(s, nchar(s), nchar(s)) == "/"){#the last character of the dir is "/"
            s = substr(s, 1, nchar(s)-1)#remove the last "/"
        }
    gdinfo = s

    ####remove those genes whose expressions are 0 across all cells
#     q = foreach(i = 1:len, .combine= c)%dopar%{
#         a = scExp[[i]] 
#         f1 = unique(a@i)+1#0-based index 
#         f2 = setdiff(1:D, f1)
#         return(f2)
#     }

#     q = foreach(a = scExp, .combine= c)%dopar%{
# #         a = scExp[[i]] 
#         f1 = unique(a@i)+1#0-based index 
#         f2 = setdiff(1:D, f1)
#         return(f2)
#     }
#     
#     f = table(q)#statistics of repeated numbers
# 
#     x = which(f==len)#all zeros in different lists
# 
#     nx = as.numeric(names(x))#selected all-zero elements
# 
#     ny = setdiff(1:D, nx)#non-zero genes
#     
#     D = length(ny)
# #     
# #     nw = list()
#     nw = foreach(i = 1:len)%dopar%{
#         z = scExp[[i]][ny, ]
#         return(z)
#     }
#     rm("scExp")
#     gc()
#     scExp = nw
#     rm("nw")
#     gc()
    
    
    
#     my = scExp
    
#     
# #     my = my[apply(my,1,function(x) sd(x)!=0),]
# #     my <- t(scale(t(my)))
# 
# #     mat = as.matrix(my)
# #     mat = my
# 
# #     D = nrow(mat)#number of genes
    
#     D = 10
#     cat("Number of valid genes:", D, "\n")
#     gnames = rownames(scExp[[1]])[ny]

    uy = y$unique_pred_clusters#unique
    if(max(as.numeric(uy)) > length(uy)){#in case we use a larger value than the number of unique elements
        y$pred_clusters  = match(y$pred_clusters, uy)#we use the index instead
    }

#     dt <- data.frame(yy = t(mat), ig = y$pred_clusters)#gene expression; cluster
    label = as.numeric(y$pred_clusters)
    rm("y")
    gc()
    dt <- data.frame(ig = label)#gene expression; cluster
   
#     outdir = "tmp"
#     if(!dir.exists(outdir)){#the folder to store the results
# 	  dir.create(outdir)
#     }else{
#         unlink(paste0(outdir, "/*"))
#     }
#     
#     pn = 1
#     k = ceiling(D/pn)
# #     D = 1000
#     foreach(i = 1:D)%dopar%{
# #          if(i!=k){
# #             gi = ((i-1)*pn+1):i*pn
# #          }else{
# #             gi = ((i-1)*pn+1):D
# #          }
#          
#          r = unlist(sapply(1:len, function(x) scExp[[x]][i,]))#the rank of the expression for the i-th gene across single cells
# #          s1 = which(r!=0)
# #          s2 = r[s1]
# #          names(s2) = s1
#          fn = paste(outdir, "/t", i, ".rds", sep = "")
#          saveRDS(r, fn)
#     }
# #     rm(scExp)
# #     rm(y)
#     gc()

#     allfiles = list.files(scExp_dir, full.names = TRUE)#include the full path
#     len = length(allfiles)
    
    allfiles = list.files(gdinfo, full.names = TRUE)#include the full path
    xy1 = as.numeric(gsub("\\D*([0-9]+).*$", "\\1", allfiles))#
    allfiles = allfiles[order(xy1)]#reorganize the file order according to the file name with digits
    glen = length(allfiles)#number of partitions
    
    cat("Checking genes...\n")

#     foreach(i = 1:len)%dopar%{
#         k = paste0("d",i)
#         assign(k, scExp[[i]])
#         
#         return(NULL)
#     }
#     rm("scExp")
#     gc()
#         
#     
#         foreach(j = 1:D)%:%foreach(i = 1:len)%dopar%{
#             k = paste0("d",i,"g",j)
#             assign(k, scExp[[i]][j,])
#             return(NULL)
#         }
#         rm("scExp")
#         gc()
    

#     D = 10
#     total <- D
#     # create progress bar
#     pb <- txtProgressBar(min = 0, max = total, style = 3)
    ######parallel
#     theta = 1e-5

#     myenv = new.env()
#     myenv = foreach(i = 1:20)%dopar%{
# #         dn = paste0("d", i)
#         unlist(sapply(1:len, function(x) scExp[[x]][i,]))
#     }
#     ginfo = matrix(0, ncol = 4, nrow = D)
#     ginfo = foreach(i = 1:D)%dopar%{
#     g5 = foreach(i = 1:len, .combine = "rbind")%dopar%{
    ginfo = foreach(i = 1:glen, .combine = "rbind")%do%{
#         tt = character(5)
        ddk = readRDS(allfiles[i])
        ng = nrow(ddk)
#     ginfo = foreach(i = 1:length(ny), .combine = c)%dopar%{
#         cat("\nGene", ny[i], "\n")
#         r = rank(mat[i, ])#gene rank
#         cat("Collect the", i, "-th gene expression across 1.3 million cells...\n")
#         eval(parse(text = paste0("mat <- myenv[[", i,"]]")))
#         fn = paste(outdir, "/t", i, ".rds", sep = "")
#         mat = readRDS(allfiles[i])
#         r = rank(mat)
#         r = rank(unlist(sapply(1:len, function(x) scExp[[x]][i,])))#the rank of the expression for the i-th gene across single cells
        
        gn = rownames(ddk)#gene name
#         tt = matrix(character(0), ncol = 4, nrow = ng)
        tt = foreach(j = 1:ng, .combine = "rbind")%dopar%{
            dd = ddk[j, ]
            dp = length(which(dd!=0))/length(dd)
#             gn = rownames(ddk)[j]#gene name
            cat("Number", i, "Gene", j, ":", gn[j], "---- Non-zero percentage:", dp, "\n")
       
            if(dp > theta){
            r = rank(dd)
#             r = dd
#         s1 = which(dd!=0)
#         dt1 = dt[s1]
#         s2 = aggregate(r1~dt1, FUN = mean)
        
#         cat("Get average expression on the", i,"-th across different clusters...\n")
   
            s = aggregate(r~dt$ig, FUN = mean)#average over one gene across cells
#         cat("Get the cell index of the maximum expression...\n")
            icluster = which.max(s$r)#the cluster with the maximum value
            x1 = r[as.numeric(dt$ig) == icluster]
            x2 = r[as.numeric(dt$ig) != icluster]
            p = wilcox.test(x1, x2)$p.value
    
            pred <- prediction(r, as.numeric(dt$ig) == icluster)
            auc0 <- unlist(performance(pred, "auc")@y.values)
        }else{
            auc0 = icluster = 0
            p = 1
        }
    
            
            
#             dmat = numeric(0)
#             if(p < pvalue && auc0 > auc){
#                 dmat = rbind(dmat, ddk[j, , drop = FALSE])
#                 f = 1
#             }else{
#                 f = 0
#                 dmat = NULL
#             }
#             
            t0 = c(auc0, icluster, p, dp)#
            rm(ddk)
            gc()
            
#             t1 = list()
#             t1$t0 = t0
#             t1$dmat = dmat
#             
            return(t0)
    
        }
        rownames(tt) = gn
#         dd = unlist(sapply(1:len, function(x) scExp[[x]][ny[i],]))
#         dp = length(which(dd!=0))/length(dd)
#         gn = rownames(scExp[[1]])[ny[i]]#gene name
#         cat("Number", i, "Gene", ny[i], "---- Non-zero percentage:", dp, "\n")
#         
        
        # update progress bar
#         setTxtProgressBar(pb, i)
    
        # cat(c(auc, icluster, p), "\n")
        return(tt)
    }
    

#     ginfo = matrix(ginfo, nrow = 3)
    ###############


#     close(pb)
    
#     ginfo = g5#collect those info
#     storage.mode(ginfo) = "numeric"#convert to numeric
    ginfo = data.frame(ginfo)
    
#     rownames(ginfo) = g5[, 5]

#     ginfo <- data.frame(matrix(ginfo, ncol = 5, byrow = T))#it is not necessarily to be a data frame
    colnames(ginfo) = c("auc", "icluster", "pvalue", "sparsity")
    
#     cat(ginfo, "\n")
    nr = nrow(ginfo)
    cat("Number of genes checked:", nr, "\n")
#     exmat = scExp[[1]]
#     rownames(ginfo) = gnames[1:nr]
    
    ginfo = ginfo[ginfo$sparsity > theta, ]#remove low-expressed genes
    
    ginfo = ginfo[!is.nan(ginfo$pvalue), ]#remove those corresponding to p-value == NaN

    ginfo$pvalue = p.adjust(ginfo$pvalue, method = "holm")

    #after p-value adjustment
    gallinfo = ginfo#note that nrow(gallinfo) may not be the number of genes, because we have removed those sparse ones
    colnames(gallinfo) = c("auc", "icluster", "pvalue", "sparsity")
    ###in case some clusters do not have any qualified marker genes, we need to lower the AUC threshold (0.85)
    dt = data.table(ginfo)
    dt0 = dt[, list(maxauc = max(auc), minpvalue = min(pvalue)), by = icluster]

    ##adjusted AUC and p-value
#     adpvalue = max(0.01, max(dt0$minpvalue))
    adpvalue = pvalue
    adauc = min(auc, min(dt0$maxauc))
#     adauc = auc

    x1 = which(ginfo$pvalue < adpvalue)

    x2 = which(ginfo$auc > adauc)

    # x1 = which(ginfo$pvalue < 0.01)
    # 
    # x2 = which(ginfo$auc > 0.85)

    xx = intersect(x1, x2)

    mginfo = ginfo[xx, ]#selected genes with p-values and their clusters
    
    nsgene = nrow(mginfo)
    cat("Number of selected marker genes:", nsgene, "\n")
    #the corresponding expressions of marker genes
#     dd = foreach(i = 1:nsgene, .combine = "rbind")%dopar%{
#         unlist(sapply(1:len, function(x) scExp[[x]][rownames(mginfo[i, ]),]))
#     }
    
#     dd = foreach(i = 1:len, .combine = "cbind")%dopar%{
#         colnames(scExp[[i]]) = NULL#remove the cell id
#         x = scExp[[i]][rownames(mginfo), , drop = FALSE]#even for one row, it is still regarded as a matrix with the nrow == 1
#         return(x)
#     }
#     dd = Matrix

#     adpvalue = max(mginfo$pvalue)
#     adauc = min(mginfo$auc)
    
    end_time <- Sys.time()
    
    t <- difftime(end_time, start_time, units = "mins")  #difference time in minutes
    cat("Running time:", t ,"minutes\n")
    cat("-----------------------------------------------------------------------\n")
    cat("The adjusted p-avalue is", adpvalue, "and the adjusted AUROC is", adauc, "\n")
    
    res = list()
    res$mginfo = mginfo
#     res$mat = dd
    res$label = label
    res$gallinfo = gallinfo
    return(res)
}

# rankcells <- function(x){
#     n = length(x)
#     s1 = which(x!=0)#indices of non-zero cells
#     s2 = setdiff(1:n, s1)#indices of zero cells
#     ns2 = length(s2)#number of zeros
#     
#     y = numeric(n)
#     y[s2] = (1+ns2)/2#ranking for those zeros
#     y[s1] = rank(x[s1]) + ns2#ranking for those non-zeros
#     
#     return(y)
# }

# c2t <- function(w, n.cores){
#     registerDoParallel(n.cores)
#     len = length(w)
#     
#     s = list()
#     s = foreach(i = 1:len)%dopar%{
#         as(w[[i]], "dgTMatrix")
#     }
#     stopImplicitCluster()
# }
# 
# sp2 <- function(s){
#     f = cbind(s@i, s@j, s@x)
#     
#     return(f)
# }

