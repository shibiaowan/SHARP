#' Run SHARP for single-cell RNA data clustering
#'
#' SHARP: Single-cell RNA-Seq Hyper-fast and Accurate clustering via ensemble Random Projection.
#'
#' @param scExp input single-cell expression matrix, where each column represents a cell and each row represents a gene.
#' @param exp.type  the data type of single-cell expression matrix. Common types include 'count', 'UMI', 'CPM', 'TPM', 'FPRKM' and 'RPKM'. If missing, SHARP regards scExp as already normalized expression matrix.
#' @param ensize.K  number of applications of random projection for ensemble. The default value is 15.
#' @param reduced.ndim  the dimension to be reduced to. If missing, the value will be estimated by an equation associated with number of cells (see our paper and supplementary materials for details).
#' @param base.ncells   a base threshold of number of cells. The default value is 5000. When the number of cells of a dataset is smaller than this threshold, we use SHARP_small function; otherwise, we use SHARP_large.
#' @param partition.ncells  number of cells for each partition when using SHARP_large. The default value is 2000.
#' @param finalN.cluster    number of clusters for the final clustering results. The default is NULL, i.e., without giving the number of clusters, and SHARP will automatically determine the optimal number of clusters. If given, SHARP will calculate according to the given number of clusters.
#' @param enpN.cluster  number of clusters for the weighted ensemble meta-clustering only for SHARP_large. The default is NULL, i.e., without giving the number of clusters, and SHARP will automatically determine the optimal number of clusters. If given, SHARP will calculate according to the given number of clusters.
#' @param indN.cluster  number of clusters for the individual RP-based hierarchical clustering. The default is NULL, i.e., without giving the number of clusters, and SHARP will automatically determine the optimal number of clusters. If given, SHARP will calculate according to the given number of clusters.
#' @param n.cores   number of cores to be used. The default is (n-1) cores, where n is the number of cores in your local computer or server.
#' @param rN.seed   a number using which we can set seeds for SHARP to obtain reproducible results.
#'
#' @details This is the main interface for SHARP to process and analyze different kinds of single-cell RNA-Seq data. Only one parameter is manadatory, i.e., scExp, the single-cell expression matrix. In most cases, most of the parameters can be determined automatically or have been optimized, so users don't have to take efforts to try different parameters. While for some other cases where users need to change parameters, SHARP also provides various parameters, including algorithm-related parameters, hierarchical-clustering-related parameters, parallel-computing parameters and parameters to obtain reproducible results, for better optimizing the performance.
#'
#' @return a list containing the SHARP clustering results, the predicted optimal number of clusters, time SHARP consumes for clustering, some intermediate results including clustering results by each random-projection based hierarchical clustering and other related statstical information including number of cells, genes, reduced dimensions and number of applications of random projection.
#'
#' @examples
#' enresults = SHARP(scExp)
#'
#' @author Shibiao Wan <shibiao@pennmedicine.upenn.edu>, Junil Kim <junilkim@pennmedicine.upenn.edu>, Kyoung Jae Won <wonk@pennmedicine.upenn.edu>
#'
#' @import foreach
#'
#' @import parallel
#'
#' @import doParallel
#'
#' @export
SHARP <- function(scExp, exp.type, ensize.K, reduced.ndim, base.ncells, partition.ncells, 
    hmethod, finalN.cluster = NULL, enpN.cluster = NULL, indN.cluster = NULL, minN.cluster, 
    maxN.cluster, sil.thre, height.Ntimes, n.cores, rN.seed) {
    # timing
    start_time <- Sys.time()  #we exclude the time for loading the input matrix
    
    title = "scRNA-Seq Clustering"
    
    if (missing(scExp)) {
        stop("No expression data is provided!")
    }
    
    ngenes = nrow(scExp)  #number of genes
    ncells = ncol(scExp)  #number of cells
    
    if (!missing(exp.type)) {
        if (exp.type == "count" || exp.type == "UMI") {
            scExp = log10(t(t(scExp)/colSums(scExp)) * 1e+06 + 1)
        } else if (exp.type == "CPM" || exp.type == "TPM" || exp.type == "FPKM" || 
            exp.type == "RPKM") {
            scExp = log10(scExp + 1)
        }
    }
    
    if (missing(reduced.ndim)) {
        # default dimensions to be reduced
        reduced.ndim = ceiling(log2(ncells)/(0.2^2))  #reduced 100 times of dimensions; about 200-dim
    }
    
    if (missing(base.ncells)) {
        # threshold of number of cells to determine whether SHARP_small or SHARP_large
        # should be used
        base.ncells = 5000
    }
    
    if (missing(partition.ncells)) {
        # number of cells for each partition group
        partition.ncells = 2000
    }
    
    if (missing(hmethod)) {
        # hierarchical clustering aggolomeration method
        hmethod = "ward.D"
    }
    
    # if no minimum number of clusters is provided
    if (missing(minN.cluster)) {
        minN.cluster = 2  #by default, we try the minimum number of clusters starting from 2
    }
    
    # if no maximum number of clusters is provided
    if (missing(maxN.cluster)) {
        maxN.cluster = 40  #by default, we try the maximum number of clusters as large as 40 or the number of cells minus 1, whichever is smaller.
    }
    
    # if no threshold for the maximum Silhouette index is provided
    if (missing(sil.thre)) {
        sil.thre = 0.35  #by default, we use 0.35 to determine whether we use Silhouette index as the criteria to determine the optimal number of clusters
    }
    
    # if no threshold for the height difference is provided
    if (missing(height.Ntimes)) {
        height.Ntimes = 2  #by default, we select the first height which is (height.Ntimes) times larger than the immediate consecutive height
    }
    
    if (missing(n.cores)) {
        # number of cores to be used, the default is to use all but one cores
        n.cores = detectCores() - 1
    }
    # registerDoMC(n.cores)
    registerDoParallel(n.cores)
    
    
    if (!missing(rN.seed)) {
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
    
    colorL <<- c("red", "purple", "blue", "yellow", "green", "orange", "brown", "gray", 
        "black", "coral", "beige", "cyan", "turquoise", "pink", "khaki", "magenta", 
        "violet", "salmon", "goldenrod", "orchid", "seagreen", "slategray", "darkred", 
        "darkblue", "darkcyan", "darkgreen", "darkgray", "darkkhaki", "darkorange", 
        "darkmagenta", "darkviolet", "darkturquoise", "darksalmon", "darkgoldenrod", 
        "darkorchid", "darkseagreen", "darkslategray", "deeppink", "lightcoral", 
        "lightcyan")
    
    
    # print(paste('For Dataset: ', dir1, sep = ''), quote = FALSE)
    
    # tcfile = 'id2celltype.txt' file1 = paste(dir1, '/', tcfile , sep = '') gtc =
    # read.delim(file1, check.names = F)#read the file containing the ground-truth
    # clusters
    
    
    
    # print(paste('Number of cells: ', ncells, sep = '')) print(paste('Number of
    # genes: ', ngenes, sep = '')) # print(paste('Ground-truth Number of clusters: ',
    # length(unique(gtc$cellType)), sep = '')) print(paste('The dimension has been
    # reduced from ', ngenes, ' to ', reduced.ndim, sep=''))
    
    cat("Number of cells: ", ncells, "\n")
    cat("Number of genes: ", ngenes, "\n")
    cat("The dimension has been reduced from ", ngenes, " to ", reduced.ndim, "\n")
    if (ncells < base.ncells) {
        # print('Using SHARP_small...')
        cat("Using SHARP_small...\n")
        if (missing(ensize.K)) {
            # default times of random projection
            ensize.K = 15  #K times of random projection
        }
        enresults = SHARP_small(scExp, ncells, ensize.K, reduced.ndim, hmethod, finalN.cluster, 
            indN.cluster, minN.cluster, maxN.cluster, sil.thre, height.Ntimes, rN.seed)
    } else {
        # print('Using SHARP_large...')
        cat("Using SHARP_large...\n")
        if (missing(ensize.K)) {
            # default times of random projection
            ensize.K = 5  #K times of random projection
        }
        enresults = SHARP_large(scExp, ncells, ensize.K, reduced.ndim, partition.ncells, 
            hmethod, finalN.cluster, enpN.cluster, indN.cluster, minN.cluster, maxN.cluster, 
            sil.thre, height.Ntimes, rN.seed)
    }
    
    end_time <- Sys.time()
    
    t <- difftime(end_time, start_time, units = "mins")  #difference time in minutes
    print(t)
    #################################### 
    enresults$N.cells = ncells
    enresults$N.genes = ngenes
    enresults$reduced.dim = reduced.ndim
    enresults$ensize.K = ensize.K
    enresults$time = t
    
    return(enresults)
}

#' Run SHARP for small-size (< 5000) single-cell RNA datasets
#'
#' For small-size (< 5000) datasets, we don't need to partition the datasets into several groups. Instead, we simply use ensemble random projection and weighted ensemble meta-clustering algorithms.
#'
#' @param scExp input single-cell expression matrix
#' @param ncells number of single cells
#' @param ensize.K number of applications of random projection for ensemble
#' @param reduced.ndim the dimension to be reduced to
#'
#' @examples
#' enresults = SHARP_small(scExp, ncells, ensize.K, reduced.ndim)
#'
#' @import foreach
#'
#' @import parallel
#'
#' @import doParallel
#'
#' @export
SHARP_small <- function(scExp, ncells, ensize.K, reduced.ndim, hmethod, finalN.cluster, 
    indN.cluster, minN.cluster, maxN.cluster, sil.thre, height.Ntimes, rN.seed) {
    enresults = list()
    
    p = reduced.ndim
    allrpinfo = foreach(k = 1:ensize.K) %dopar% {
        # print(paste('Random Projection: ', k, sep = ''))
        cat("Random Projection: ", k, "\n")
        # print(paste('The ', k, '-th time of random projection', sep=''), quote = FALSE)
        scExp = data.matrix(scExp)  #convert from data frame to normal matrix
        if (rN.seed == 0.5) {
            # print('working!')
            newE = RPmat(scExp, p, 0.5)  #do the RP;the result is a list#for reproducible results
        } else {
            newE = RPmat(scExp, p, 50 + rN.seed + k)  #do the RP;the result is a list#for random results
        }
        # newE = RPmat(scExp, p, seed.logic)#do the RP;the result is a list
        E1 = t(newE$projmat)  #for those which need tranpose
        
        tag = paste("_RP", p, "_", k, sep = "")  #k is the application times of random projection; p is the reduced dimension
        tmp = getrowColor(E1, hmethod, indN.cluster, minN.cluster, maxN.cluster, 
            sil.thre, height.Ntimes)
        rowColor = tmp$rowColor
        # rowColor= getrowColor(E1, tag, outdir, colorL)#hierarchical clustering metrics=
        # ARI(gtc, rowColor)#performance evaluation print(metrics)
        
        # enrp[,k] = rowColor
        
        rpinfo = list()
        # for different parameters, we do not need to save individual randome matrices
        # rpinfo$rpmat = E1#the after-random-projected matrix rpinfo$R = newE$R#the
        # random matrix
        rpinfo$tag = tag  #tag
        rpinfo$rowColor = rowColor  #the resulting clusters
        # rpinfo$metrics = metrics#the performance for each individual RPs
        rpinfo$N.cluster = length(unique(rowColor))
        return(rpinfo)
        # rpname = paste('RP_', k, sep = '') allrpinfo[[rpname]] = rpinfo
    }
    
    z = length(allrpinfo)
    enrp = matrix("0", nrow = ncells, ncol = ensize.K)
    enE <- matrix(0, nrow = ncells, ncol = p)  #already tranposed; reshuffled matrix
    for (j in 1:z) {
        # partition
        z1 = allrpinfo[[j]]
        enrp[, j] = z1$rowColor
    }
    
    finalrowColor = wMetaC(enrp, hmethod, enN.cluster = finalN.cluster, minN.cluster, 
        maxN.cluster, sil.thre, height.Ntimes)
    # finalmetrics = ARI(gtc, finalrowColor)#performance evaluation print('The
    # ensemble performance metrics are:', quote = FALSE) print(finalmetrics)
    
    
    # save(enrp, file=paste(outdir,'enrp_', K, 'times.RData', sep = ''))
    # save(finalrowColor, file = paste(outdir,'finalrowColor_', K, 'times.RData', sep
    # = '')) save(finalmetrics, file = paste(outdir,'finalmetrics_', K,
    # 'times.RData', sep = ''))
    
    # enresults$enrp = enrp#we have already saved the rowColor for each individual
    # RP, thus we do not need to save enrp
    
    # enresults$finalrowColor = finalrowColor$finalC
    enresults$final_pred_clusters = finalrowColor$finalC
    enresults$unique_pred_clusters = unique(finalrowColor$finalC)
    # enresults$finalmetrics = finalmetrics
    enresults$N.pred_cluster = length(unique(finalrowColor))
    enresults$allrpinfo = allrpinfo
    # save(enresults, file=paste(outdir,'enresults_', p, 'dim_', K, 'times.RData',
    # sep = ''))
    return(enresults)
}

#' Run SHARP for large-size (>= 5000) single-cell RNA datasets
#'
#' For large-size (>= 5000) datasets, we suggest first partitioning the datasets into several groups, then we run SHARP for each group, and finally and we ensemble the results of each group by a similarity-based meta-clustering algorithm.
#'
#' For each partition (or group), the default number of cells is set to 2000 for each group. The users can also set a different number according to the computational capability of their own local computers. The suggested criteria to set this number is that as long as SHARP_small can run in a fast enough (depending on users' requirements) way for the selected number of single cells.
#'
#' @param scExp input single-cell expression matrix
#' @param ncells number of single cells
#' @param ensize.K number of applications of random projection for ensemble
#' @param reduced.dim the dimension to be reduced to
#' @param partition.ncells number of cells for each partition when using SHARP_large
#'
#' @examples
#' enresults = SHARP_large(scExp, ncells, ensize.K, reduced.dim, partition.ncells)
#'
#' @import foreach
#'
#' @import parallel
#'
#' @import doParallel
#'
#' @export
SHARP_large <- function(scExp, ncells, ensize.K, reduced.dim, partition.ncells, hmethod, 
    finalN.cluster, enpN.cluster, indN.cluster, minN.cluster, maxN.cluster, sil.thre, 
    height.Ntimes, rN.seed) {
    # print('The Divide-and-Conquer Strategy is selected!')
    cat("The Divide-and-Conquer Strategy is selected!\n")
    ######## Partition the large data into several groups###########
    p = reduced.dim
    entag = paste("_enRP", p)
    
    K = ensize.K
    allrpinfo <- vector("list", length = K)  #declare a matrix of lists
    enresults = list()
    # colnames(rpinfo) = c('E', 'tag', 'rowColor', 'metrics')
    enrp <- matrix("0", ncells, K)  #the ensemble results after several random projection;namely several rowColor's
    
    if (rN.seed == 0.5) {
        # for reproducible results
        reind = sample(ncells)  #randomly reshuffle the data
    } else {
        set.seed(50)
        reind = sample(ncells)  #randomly reshuffle the data, but reproducible
    }
    
    # reshuffle the data
    E = scExp
    E = E[, reind]
    # ng = 2000#number of cells for each group
    
    # if(ncells < 10000){#if it is a small dataset, simply divide it into 5 parts ng
    # = ceiling(ncells/5) }
    ng = partition.ncells
    T = ceiling(ncells/ng)  #number of partitions/groups
    
    # resT = T*ng - nrow(E1)#we need to fill resT 0 for consistency and ease of
    # following processing
    if (T > 1) {
        folds <- cut(seq(1, T * ng), breaks = T, labels = FALSE)  #note that the index may be larger than the number of cells
        
        # if(ncells - (T-1)*ng < 41){
        nt = ncells - (T - 2) * ng
        nind = which(folds == (T - 1))
        folds[nind[floor(nt/2) + 1:ng]] = T  #To avoid imbalanced problems, the last two folds are roughly averaged.
        folds = folds[1:ncells]
        # } print(paste('The number of cells in Fold ', names(table(folds)), 'are: ',
        # table(folds)))
        cat("The number of cells in Folds ", names(table(folds)), "are: ", table(folds), 
            "\n")
        
        np = table(folds)
    } else {
        folds = rep(1, each = ng)
        folds = folds[1:ncells]
        np = length(folds)
    }
    
    # preparing K RP matrices
    rM = foreach(k = 1:K) %dopar% {
        if (rN.seed == 0.5) {
            ranM(E, p, 0.5)  #for random results
        } else {
            ranM(E, p, 50 + rN.seed + k)  #for reproducible results
        }
        
    }
    
    #### parallel programming#### rerowColor = foreach(t=1:T, .combine = 'c') %dopar%{
    enlist = foreach(k = 1:K, .combine = rbind) %:% foreach(t = 1:T) %dopar% {
        # nested looping; the first for different applications of random projection; the
        # second for different partitions print(paste('The ', k, '-th random projection;
        # the ', t, '-th partition', sep=''))
        tind = which(folds == t, arr.ind = TRUE)  #find the indices
        # if(t == T){ tind = tind[tind<=ncells]#select only those within the indices }
        
        # print(paste('Random Projection: ', k, ', Fold: ', t, ', Cell Number: ',
        # length(tind), sep = ''))
        cat("Random Projection: ", k, ", Fold: ", t, ", Cell Number: ", length(tind), 
            "\n")
        ######## Cluster each group########
        newE = E[, tind]  #the matrix for each group
        newE = log10(newE + 1)#logarithm transform
        
        inE = data.matrix(newE)
        
        # using the k-th random matrix
        E1 = 1/sqrt(p) * t(rM[[k]]) %*% inE
        # newinE = RPmat(inE, p)#do the RP;the result is a list E1 =
        # t(newinE$projmat)#for those which need tranpose
        
        E1 = t(E1)
        
        newE1 = data.matrix(E1)  #convert the spase dim-reduced matrix back to full dim-reduced matrix
        # if(ncells>10000){#for large-scale datasets, we normalize the values for
        # memory-efficient computation E1 = round(E1/100, digits = 3)#due to the high
        # dimensions (>10000) in the original data, the values of E1 may be very large,
        # better to normalize to a small scale for later clustering }
        
        tmp = getrowColor(newE1, hmethod, indN.cluster, minN.cluster, maxN.cluster, 
            sil.thre, height.Ntimes)  #hierarchical clustering; the result is a list containing both the predicted clusters and the maximum silhouette index
        # metrics= ARI(gtc[reind[tind], ], tmp$rowColor)#performance evaluation
        # print(metrics) rerowColor[, t] = paste(tmp, '_', t, sep = '')#for
        # distinguishing for each smaller group clustering rind = c((t-1)*ng + 1:
        # (t-1)*ng + length(tind))
        ptc = paste(tmp$rowColor, "p", t, sep = "")  #using the letter 'p' to represent the 'partial', t is the t-th part of the data
        # print(paste('ptc length is: ', length(ptc), sep = ''))
        
        tmplist = list()
        tmplist$rpc = ptc
        tmplist$pE1 = newE1
        tmplist$ind = c(k, t)
        tmplist$maxsil = tmp$maxsil
        return(tmplist)
        # rerowColor[rind] = paste(tmp, '_', t, sep = '')#for distinguishing for each
        # smaller group clustering
    }
    
    
    z = dim(enlist)
    Elist = vector("list", length = T)
    enE <- matrix(0, nrow = ncells, ncol = p)  #already tranposed; reshuffled matrix
    for (j in 1:z[2]) {
        # partition
        eind = which(folds == j)
        for (i in 1:z[1]) {
            # applications of random projection
            z1 = enlist[i, j]
            nz = names(z1)
            ptc = z1[[nz]]$rpc  #the clustering results by RP for a partitio of the data
            pE = z1[[nz]]$pE1  #the partial dim-reduced matrix
            maxsil = z1[[nz]]$maxsil
            
            
            enrp[eind, i] = ptc  #the clustering results of different random projection
            # print(paste('enrp[folds == j, i] length is: ', length(enrp[folds == j, i]), sep
            # = ''))
            
            # enE[folds == j, ] = enE[folds == j, ] + pE*maxsil
            
            enE[eind, ] = enE[eind, ] + pE
            
        }
    }
    
    
    # finalrowColor = wMetaC(enrp) parallel programming####
    frowColor = foreach(t = 1:T, .combine = "c") %dopar% {
        tind = which(folds == t, arr.ind = TRUE)
        # tind = seq((t-1)*ng + 1, t*ng) if(t == T){ tind = tind[tind<=ncells] # tind =
        # seq((T-1)*ng + 1, ncells)#select only those within the indices }
        ftmp = wMetaC(enrp[tind, ], hmethod, enpN.cluster, minN.cluster, maxN.cluster, 
            sil.thre, height.Ntimes)  #ensemble of several applications of RPs
        # enmetrics = ARI(gtc[reind[tind], ], ftmp$finalC)#performance evaluation
        # print(enmetrics) ftmp$finalC
        paste(ftmp$finalC, "en", t, sep = "")  #the letters 'en' represent the ensemble clustering
    }
    # stop('Trial 1 success!') frowColor = unlist(sapply(1:T, function(t){tind =
    # seq((t-1)*ng + 1, t*ng); if(t == T){tind = seq((T-1)*ng + 1, ncells)}; ftmp =
    # wMetaC(enrp[tind,]); paste(ftmp$finalC, '_', t, sep = '')}))
    
    # reorganizing the total cells
    if (T == 1) {
        # if the dataset is not partitioned
        SrowColor = frowColor
        N.cluster = length(unique(SrowColor))  #note that the number of clusters is determined by the unique number in the final round.
        print(paste("The optimal number of clusters is: ", N.cluster, sep = ""))
        cat("The optimal number of clusters is: ", N.cluster, "\n")
    } else {
        # newinE = RPmat(data.matrix(E), p)#do the RP;the result is a list E1 =
        # t(newinE$projmat)#for those which need tranpose E1 = data.matrix(E1)
        E1 = enE/K
        
        stmp = sMetaC(frowColor, E1, folds, hmethod, finalN.cluster, minN.cluster, 
            maxN.cluster, sil.thre, height.Ntimes)
        SrowColor = stmp$finalColor
        # SrowColor = sMetaC(frowColor, E1, folds) SrowColor = sMetaC(frowColor, E, p,
        # folds)
    }
    finalrowColor = vector(mode = "character", length = length(SrowColor))  #initialization
    finalrowColor[reind] = SrowColor  #reorganizing the final results
    # finalmetrics = ARI(gtc, finalrowColor)#performance evaluation #
    # print(paste('The predicted number of clusters is: ',
    # length(unique(finalrowColor)), sep = '')) print('The ensemble performance
    # metrics are:') print(finalmetrics)
    
    # #visualization enE = enE/K#averaging the projected matrix for visualization
    # rtsne_out <- Rtsne(as.matrix(enE), check_duplicates = FALSE) file_plot <-
    # paste('vi_',dir1, '.png', sep = '') png(file_plot, width = 900, height = 900)
    # plot(rtsne_out$Y, asp = 1, pch = 20, col = gtc$cellType, cex = 0.75, cex.axis =
    # 1.25, cex.lab = 1.25, cex.main = 1.5, xlab = 't-SNE dimension 1', ylab = 't-SNE
    # dimension 2', main = '2D t-SNE projection') dev.off()
    
    
    # # vAA = ftmp$pcaAA # t1 = factor(gtc$cellType)#rowColor is a vector like 'red
    # red purple brown ...'  # levels(t1) = c(1:length(levels(t1)))#convert the
    # categorical to numeric # t1 = as.numeric(as.character(t1)) # tcol =
    # cL[t1]#color # plot(vAA, col = tcol)#visualization file_plot <- paste('SHARP_',
    # dir1, '.png', sep = '') png(file_plot, width = 900, height = 900) # rtsne_out
    # <- Rtsne(as.matrix(ftmp$x0), pca = FALSE, verbose = TRUE, check_duplicates =
    # FALSE) rtsne_out <- Rtsne(as.matrix(ftmp$x0), pca = FALSE, check_duplicates =
    # FALSE) # plot(rtsne_out$Y, asp = 1, pch = 20, col = gtc$cellType, cex = 0.75,
    # cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, xlab = 't-SNE dimension 1',
    # ylab = 't-SNE dimension 2', main = '2D t-SNE projection') plot(rtsne_out$Y, asp
    # = 1, pch = 20, col = gtc$cellType, cex = 0.75, cex.axis = 1.25, cex.lab = 1.25,
    # cex.main = 1.5, xlab = 'SHARP Dim-1', ylab = 'SHARP Dim-2', main = 'SHARP
    # Visualization') dev.off() # save(enrp, file=paste(outdir,'enrp_', K,
    # 'times.RData', sep = ''))
    
    # save(finalrowColor, file = paste(outdir,'finalrowColor_', K, 'times.RData', sep
    # = '')) save(finalmetrics, file = paste(outdir,'finalmetrics_', K,
    # 'times.RData', sep = ''))
    
    
    # enresults$finalrowColor = finalrowColor
    enresults$final_pred_clusters = finalrowColor
    enresults$unique_pred_clusters = unique(finalrowColor)
    # enresults$finalmetrics = finalmetrics
    enresults$N.pred_cluster = length(unique(finalrowColor))
    
    # save(enresults, file=paste(outdir,'enresults_', K, 'times.RData', sep = ''))
    # saveRDS(enresults, file=paste(outdir,'largeresults_', dir1, '.rds', sep = ''))
    return(enresults)
}
