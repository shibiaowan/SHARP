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
#' @param N.cluster    number of clusters for the final clustering results. The default is NULL, i.e., without giving the number of clusters, and SHARP will automatically determine the optimal number of clusters. If given, SHARP will calculate according to the given number of clusters.
#' @param enpN.cluster  number of clusters for the weighted ensemble meta-clustering only for SHARP_large. The default is NULL, i.e., without giving the number of clusters, and SHARP will automatically determine the optimal number of clusters. If given, SHARP will calculate according to the given number of clusters.
#' @param indN.cluster  number of clusters for the individual RP-based hierarchical clustering. The default is NULL, i.e., without giving the number of clusters, and SHARP will automatically determine the optimal number of clusters. If given, SHARP will calculate according to the given number of clusters.
#' @param minN.cluster  the minimum number of clusters that SHARP will try when determining the optimal number of clusters
#' @param maxN.cluster  the maximum number of clusters that SHARP will try when determining the optimal number of clusters
#' @param sil.thre  the threshold of Silhouette index that SHARP will use the Silhouette index to determine the optimal number of clusters. In other words, if the maximum Silhouette index is larger than sil.thre, then SHARP uses the Silhouette index to determine the number of clusters; otherwise, SHARP uses the other indices (i.e., CH index and/or hierarchical heights) to determine 
#' @param height.Ntimes the number of times of the height versus the immediate next height in the hierarchical clustering. SHARP uses this parameter as a threshold to determine the location where to cut the hierarchical tree. In other words, if the current height is (height.Ntimes) times larger than the immediate next height in the descending order of heights, then SHARP cuts the tree at the median of these two heights.
#' @param logflag   a logical to determine whether to check a log-transform of the input expression matrix. By default, logflag = TRUE, i.e., SHARP will check the log-transform operation.
#' @param sncells   number of cells randomly selected for checking log-transform is necessary or not. By default, sncells = 100.
#' @param n.cores   number of cores to be used. The default is (n-1) cores, where n is the number of cores in your local computer or server.
#' @param forview   a logical to indicate whether those feature-vectors for data visualization should be saved or not. By default, it is TRUE.
#'
#' @param rM if provided, it should be a list of random matrices for random projection; otherwise, it will be calculated by SHARP_large.
#
#' @param rN.seed   a number using which we can set seeds for SHARP to obtain reproducible results.
#'
#' @details This is the main interface for SHARP to process and analyze different kinds of single-cell RNA-Seq data. Only one parameter is manadatory, i.e., scExp, the single-cell expression matrix. In most cases, most of the parameters can be determined automatically or have been optimized, so users don't have to take efforts to try different parameters. While for some other cases where users need to change parameters, SHARP also provides various parameters, including algorithm-related parameters, hierarchical-clustering-related parameters, parallel-computing parameters and parameters to obtain reproducible results, for better optimizing the performance.
#'
#' @return a list containing the SHARP clustering results, distribution of the clustering results, the predicted optimal number of clusters, time SHARP consumes for clustering, some intermediate results including clustering results by each random-projection based hierarchical clustering and other related statstical information including number of cells, genes, reduced dimensions and number of applications of random projection.
#'
#' @examples
#' enresults = SHARP(scExp)
#'
#' @author Shibiao Wan <shibiao@pennmedicine.upenn.edu>
#'
#' @import foreach
#'
#' @import parallel
#'
#' @import doParallel
#'
#' @export
SHARP <- function(scExp, exp.type, ensize.K, reduced.ndim, base.ncells, partition.ncells, 
    hmethod, N.cluster = NULL, enpN.cluster = NULL, indN.cluster = NULL, minN.cluster, 
    maxN.cluster, sil.thre, height.Ntimes, logflag, sncells, n.cores, forview = TRUE, prep, rM, rN.seed) {
    # timing
    start_time <- Sys.time()  #we exclude the time for loading the input matrix
    
    title = "scRNA-Seq Clustering"
    
    if (missing(scExp)) {#scRNA-seq expression matrix
        stop("No expression data is provided!")
    }
    
    ngenes = nrow(scExp)  #number of genes
    ncells = ncol(scExp)  #number of cells
    cat("-----------------------------------------------------------------------\n")
    cat("Data info:\n")
    cat("Number of cells:", ncells, "\n")
    cat("Number of genes:", ngenes, "\n")
    
#     if (!missing(exp.type)) {
#         if (exp.type == "count" || exp.type == "UMI") {
#             scExp = log10(t(t(scExp)/colSums(scExp)) * 1e+06 + 1)
#         } else if (exp.type == "CPM" || exp.type == "TPM" || exp.type == "FPKM" || 
#             exp.type == "RPKM") {
#             scExp = log10(scExp + 1)
#         }
#     }
    
    cat("-----------------------------------------------------------------------\n")
    cat("Preprocessing:\n")
    if(missing(prep)){
        if(ncells < 1e4){
            prep = TRUE
        }else{
            prep = FALSE
        }
    }
    
    if(prep){
        if(any(scExp < 0)){
            warning("Your expression matrix contain negative values! SHARP will replace negative values with 0!\n")
            scExp[scExp < 0] = 0
        }
    
        scExp = scExp[rowSums(scExp) != 0, ]#remove those all-zero genes
    }
    
    
    cat("Normalization...\n")
    if (!missing(exp.type)) {#expression type
        cat("Converting expression type...\n")
        if (exp.type != "CPM" && exp.type != "TPM") {#if the type is neither CPM nor TPM
            scExp = t(t(scExp)/colSums(scExp)) * 1e+06#we do normalization
        }
    }else{
        # scExp = t(t(scExp)/colSums(scExp)) * 1e+06#we do normalization
    }
   
    if (missing(reduced.ndim)) {#reduced dimension
        # default dimensions to be reduced
        reduced.ndim = ceiling(log2(ncells)/(0.2^2))  #reduced 100 times of dimensions; about 200-dim
    }
    
    if (missing(base.ncells)) {
        # threshold of number of cells to determine whether SHARP_small or SHARP_large
        # should be used
        base.ncells = 5000
    }
    
    if (missing(partition.ncells)) {# number of cells for each partition group
        partition.ncells = 2000
    }
    
    if (missing(hmethod)) {# hierarchical clustering aggolomeration method
        hmethod = "ward.D"
    }
    
    # if no minimum number of clusters is provided
    if (missing(minN.cluster)) {
        minN.cluster = 2  #by default, we try the minimum number of clusters starting from 2
    }
    
    # if no maximum number of clusters is provided
    if (missing(maxN.cluster)) {
        maxN.cluster = max(40, ceiling(ncells/5000))  #by default, we try the maximum number of clusters as large as 40 or the number of cells minus 1, whichever is smaller.
    }
    
    # if no threshold for the maximum Silhouette index is provided
    if (missing(sil.thre)) {
        sil.thre = 0.35  #by default, we use 0.35 to determine whether we use Silhouette index as the criteria to determine the optimal number of clusters
    }
    
    # if no threshold for the height difference is provided
    if (missing(height.Ntimes)) {
        height.Ntimes = 2  #by default, we select the first height which is (height.Ntimes) times larger than the immediate consecutive height
    }
    
    if(missing(rM)){
        rM = TRUE
    }
    
    if (missing(n.cores)) {
        # number of cores to be used, the default is to use all but one cores
        n.cores = detectCores() - 1
    }
    # registerDoMC(n.cores)
#     registerDoParallel(n.cores)
    
    if (!missing(rN.seed)) {#random seed for reproducible results
        # seed number for reproducible results
        if (!is.numeric(rN.seed)) {
            stop("The rN.seed should be a numeric!")
        }
        if (rN.seed%%1 != 0 && rN.seed != 0.5) {
            stop("The rN.seed should be an integer!")
        }
    } else {
        rN.seed = 0.5  #use a non-integer to represent random results
    }
    
    if (!missing(N.cluster)){#if the number of clusters is pre-defined, we adopted divide-and-conquer strategies regardless of the number of single cells
#         enpN.cluster = N.cluster#for ensemble clustering
        if (ncells < base.ncells){
            indN.cluster = N.cluster
            base.ncells= ceiling(ncells/2)
            partition.ncells=ceiling(ncells/2)
            if (missing(ensize.K)) {
                ensize.K = 15
            }
        }
    }
    

    
#     colorL = c("purple",  "pink",  "black",  "orange", "turquoise", "yellow", "beige", "gray", 
#          "coral",    "khaki",  "violet", "magenta",
#          "salmon", "goldenrod", "orchid", "seagreen", "slategray", "darkred", 
#         "darkblue", "darkcyan", "darkgreen", "darkgray", "darkkhaki", "darkorange", 
#         "darkmagenta", "darkviolet", "darkturquoise", "darksalmon", "darkgoldenrod", 
#         "darkorchid", "darkseagreen", "darkslategray", "deeppink", "lightcoral", 
#         "lightcyan")
    if (missing(logflag)){#whether we need to check log-transform; 
        if(ncells < 1e4){
            logflag = TRUE
        }else{#if the number of cells is too large, we do not check but directly adopt the log-transform
            logflag = FALSE
        }
        
    }
    
    if (logflag){#if we want to check whether log-transform or not
        if (missing(sncells)){#number of cells selected for determining whether log transform or not
            #by default, we randomly select 100 cells to determine log transform or not
            sncells = 100
        }
        flag = testlog(scExp, ncells, reduced.ndim, sncells, n.cores)#test whether log transform is necessary or not; by default, it will randomly select 100 cells to test
        
        # better not directly do this, because it may involve huge computations for log transform
        if (flag){#if log transform is better
            # scExp = log2(scExp + 1)#better not directly do this, because it may involve huge computations for log transform; do it later
            cat("Log-transform is necessary!\n")
        }else{
            cat("Log-transform is not necessary!\n")
        }
    }else{
        cat("Log-transform is employed!\n")
        flag = TRUE
    }
   
    # print(paste('For Dataset: ', dir1, sep = ''), quote = FALSE)
    
    # tcfile = 'id2celltype.txt' file1 = paste(dir1, '/', tcfile , sep = '') gtc =
    # read.delim(file1, check.names = F)#read the file containing the ground-truth
    # clusters
    
    
    
    # print(paste('Number of cells: ', ncells, sep = '')) print(paste('Number of
    # genes: ', ngenes, sep = '')) # print(paste('Ground-truth Number of clusters: ',
    # length(unique(gtc$cellType)), sep = '')) print(paste('The dimension has been
    # reduced from ', ngenes, ' to ', reduced.ndim, sep=''))
    cat("-----------------------------------------------------------------------\n")
    cat("Parameter Setting:\n")
    # cat("The dimension has been reduced from ", ngenes, " to ", reduced.ndim, "\n")
    if (ncells < base.ncells) {
        # print('Using SHARP_small...')
        cat("Using SHARP_small...\n")
        if (missing(ensize.K)) {
            # default times of random projection
            ensize.K = 15  #K times of random projection
        }
        cat("Ensemble size:", ensize.K, "\n")
        cat("The dimension has been reduced from", ngenes, "to", reduced.ndim, "\n")
        cat("No partition is required!\n")
        cat("-----------------------------------------------------------------------\n")
        cat("Analysis starts...\n")
        enresults = SHARP_small(scExp, ncells, ensize.K, reduced.ndim, hmethod, N.cluster, 
            indN.cluster, minN.cluster, maxN.cluster, sil.thre, height.Ntimes, flag, n.cores, forview, rN.seed)
    } else {
        # print('Using SHARP_large...')
        cat("Using SHARP_large...\n")
        if (missing(ensize.K)) {
            # default times of random projection
            ensize.K = 5  #K times of random projection
        }
        cat("Ensemble size:", ensize.K, "\n")
        cat("The dimension has been reduced from", ngenes, "to", reduced.ndim, "\n")
        cat("Partition block size:", partition.ncells, "\n")
        cat("-----------------------------------------------------------------------\n")
        cat("Analysis starts...\n")
        enresults = SHARP_large(scExp, ncells, ensize.K, reduced.ndim, partition.ncells, 
            hmethod, N.cluster, enpN.cluster, indN.cluster, minN.cluster, maxN.cluster, 
            sil.thre, height.Ntimes, flag, n.cores, forview, rM, rN.seed)
    }
    
    end_time <- Sys.time()
    
    t <- difftime(end_time, start_time, units = "mins")  #difference time in minutes
    # print(t)
    
    cat("Analysis complete!\n")
    cat("-----------------------------------------------------------------------\n")
    cat("Total running time:", t, "minutes\n")
    #################################### 
    enresults$N.cells = ncells
    enresults$N.genes = ngenes
    enresults$reduced.dim = reduced.ndim
    enresults$ensize.K = ensize.K
    enresults$time = t
    
    ###parameter saving###
    paras = list()
    paras$ensize.K = ensize.K
    paras$reduced.ndim = reduced.ndim
    paras$base.ncells = base.ncells
    paras$partition.ncells = partition.ncells
    paras$hmethod = hmethod
    paras$N.cluster = N.cluster
    paras$minN.cluster = minN.cluster
    paras$maxN.cluster = maxN.cluster
    paras$sil.thre = sil.thre
    paras$height.Ntimes = height.Ntimes
    paras$n.cores = n.cores
    
    enresults$paras = paras
    
#     stopImplicitCluster()
    
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
SHARP_small <- function(scExp, ncells, ensize.K, reduced.ndim, hmethod, N.cluster, 
    indN.cluster, minN.cluster, maxN.cluster, sil.thre, height.Ntimes, flag, n.cores, forview, rN.seed) {
    enresults = list()
    
    if(flag){
        scExp = log10(scExp + 1)#logarithm transform
    }
    
    K = ensize.K
    p = reduced.ndim
    registerDoParallel(n.cores)
    allrpinfo = foreach(k = 1:ensize.K) %dopar% {
#         print(paste('Random Projection: ', k, sep = ''))
        pid = Sys.getpid()
        cat("Process ", pid, "----Random Projection: ", k, " out of ", K, "\n", sep = "")
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
        rpinfo$indE = E1
        return(rpinfo)
        # rpname = paste('RP_', k, sep = '') allrpinfo[[rpname]] = rpinfo
    }
    stopImplicitCluster()
    
    z = length(allrpinfo)
    enrp = matrix("0", nrow = ncells, ncol = ensize.K)
    enE <- matrix(0, nrow = ncells, ncol = p)  #already tranposed; reshuffled matrix
    for (j in 1:z) {
        # partition
        z1 = allrpinfo[[j]]
        enrp[, j] = z1$rowColor
        
        enE = enE + z1$indE
    }
    
    finalrowColor = wMetaC(enrp, hmethod, enN.cluster = N.cluster, minN.cluster, 
        maxN.cluster, sil.thre, height.Ntimes)
    # finalmetrics = ARI(gtc, finalrowColor)#performance evaluation print('The
    # ensemble performance metrics are:', quote = FALSE) print(finalmetrics)
    
    
    # save(enrp, file=paste(outdir,'enrp_', K, 'times.RData', sep = ''))
    # save(finalrowColor, file = paste(outdir,'finalrowColor_', K, 'times.RData', sep
    # = '')) save(finalmetrics, file = paste(outdir,'finalmetrics_', K,
    # 'times.RData', sep = ''))
    
    # enresults$enrp = enrp#we have already saved the rowColor for each individual
    # RP, thus we do not need to save enrp
    
    viE = as.matrix(enE/K)
    ####reorganize the final prediction results
    finalrowColor = as.numeric(finalrowColor)
    y = finalrowColor$finalC
    uy = unique(y)
    newy = match(y, uy)#we use the index instead
    newuy = unique(newy)
    
    tn = table(newy)
    tn1 = tn[order(as.numeric(names(tn)))]#order the distributions of the clustering results
    
    # enresults$finalrowColor = finalrowColor$finalC
    enresults$pred_clusters = newy
    enresults$unique_pred_clusters = sort(as.numeric(newuy))
    enresults$distr_pred_clusters = tn1
    # enresults$finalmetrics = finalmetrics
    enresults$N.pred_cluster = length(newuy)
    enresults$allrpinfo = allrpinfo
    if(forview){#if we want to visualiza data, we need to save the dimension-reduced matrices/feature vectors
        enresults$x0 = finalrowColor$x0
        enresults$viE = viE
    }
    
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
    N.cluster, enpN.cluster, indN.cluster, minN.cluster, maxN.cluster, sil.thre, 
    height.Ntimes, flag, n.cores, forview, rM, rN.seed) {
    # print('The Divide-and-Conquer Strategy is selected!')
    # cat("The Divide-and-Conquer Strategy is selected!\n")
    ######## Partition the large data into several groups###########
    p = reduced.dim
    entag = paste("_enRP", p)
    
    K = ensize.K
    allrpinfo <- vector("list", length = K)  #declare a matrix of lists
   
    # colnames(rpinfo) = c('E', 'tag', 'rowColor', 'metrics')
#     enrp <- matrix("0", nrow = ncells, ncol = K)  #the ensemble results after several random projection;namely several rowColor's
    
    if (rN.seed == 0.5) {
        # for reproducible results
        reind = sample(ncells)  #randomly reshuffle the data
    } else {
        set.seed(50)
        reind = sample(ncells)  #randomly reshuffle the data, but reproducible
    }
    
    # reshuffle the data
    E = scExp
    rm(scExp)#remove unnecessary objects for saving memory
    if(ncol(E) < 1e5){
        cat("Reshuffling the order of single cells...\n")
        E = E[, reind]
    }#if the number is larger than 50000, we do not need to reshuffle
    
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
        cat("The number of cells in Folds", names(table(folds)), "are:", table(folds), 
            "\n")
        
        np = table(folds)
    } else {
        folds = rep(1, each = ng)
        folds = folds[1:ncells]
        np = length(folds)
    }
    
    # preparing K RP matrices
    if(typeof(rM) != "list"){
        registerDoParallel(n.cores)
        rM = foreach(k = 1:K) %dopar% {
            if (rN.seed == 0.5) {
                ranM(E, p, 0.5)  #for random results
            } else {
                ranM(E, p, 50 + rN.seed + k)  #for reproducible result
            }
        }
        stopImplicitCluster()
    }
    
    
    #### parallel programming#### rerowColor = foreach(t=1:T, .combine = 'c') %dopar%{
    registerDoParallel(n.cores)
    enlist = foreach(k = 1:K, .combine= "c") %:% foreach(t = 1:T) %dopar% {
        # nested looping; the first for different applications of random projection; the
        # second for different partitions print(paste('The ', k, '-th random projection;
        # the ', t, '-th partition', sep=''))
        tind = which(folds == t, arr.ind = TRUE)  #find the indices
        # if(t == T){ tind = tind[tind<=ncells]#select only those within the indices }
        
        # print(paste('Random Projection: ', k, ', Fold: ', t, ', Cell Number: ',
        # length(tind), sep = ''))
        pid = Sys.getpid()
        cat("Process ", pid, "----Random Projection: ", k, " out of ", K, ", Fold: ", t, " out of ", T, ", Cell Number: ", length(tind), 
            "\n", sep = "")
        ######## Cluster each group########
        newE = E[, tind]  #the matrix for each group
#         cat("newE completed for:", k, "-th RP and", t, "-th partition\n")
        if(flag){#whether log-transform
            newE = log10(newE + 1)#logarithm transform
        }
# # #         newE = log10(newE + 1)#logarithm transform
        
        inE = data.matrix(newE)
        
#         cat("inE completed for:", k, "-th RP and", t, "-th partition\n")
        
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
        
#         cat("Dim reduction completed for:", k, "-th RP and", t, "-th partition\n")
        tmp = getrowColor(newE1, hmethod, indN.cluster, minN.cluster, maxN.cluster, 
            sil.thre, height.Ntimes)  #hierarchical clustering; the result is a list containing both the predicted clusters and the maximum silhouette index
        # metrics= ARI(gtc[reind[tind], ], tmp$rowColor)#performance evaluation
        # print(metrics) rerowColor[, t] = paste(tmp, '_', t, sep = '')#for
        # distinguishing for each smaller group clustering rind = c((t-1)*ng + 1:
        # (t-1)*ng + length(tind))
#         cat("hierarchical clusteirng completed for:", k, "-th RP and", t, "-th partition\n")
        ptc = paste(tmp$rowColor, "p", t, sep = "")  #using the letter 'p' to represent the 'partial', t is the t-th part of the data
        
#         print(paste('ptc length is: ', length(ptc), " for: ", k, "-th RP and ", t, "-th partition", sep = ''))
        
        tmplist = list()
        tmplist$rpc = ptc
        tmplist$pE1 = newE1
        tmplist$ind = c(k, t)
        tmplist$maxsil = tmp$maxsil
        
#         cat(str(tmplist), "for:", k, "-th RP and", t, "-th partition\n")
        
#         if(length(tmplist) == 0){
#             cat("tmplist is null!\n")
#         }
        return(tmplist)
        # rerowColor[rind] = paste(tmp, '_', t, sep = '')#for distinguishing for each
        # smaller group clustering
    }
    stopImplicitCluster()
    
#     z = dim(enlist)
    z = length(enlist)
#     cat("The length of enlist is", z, "\n")
    
#     cat("The dimension of z is", z, "\n")
#     Elist = vector("list", length = T)
    enrp <- matrix("0", nrow = ncells, ncol = K)  #the ensemble results after several random projection;namely several rowColor's
    enE <- matrix(0, nrow = ncells, ncol = p)  #already tranposed; reshuffled matrix
    for(i in 1:z){
        z1 = enlist[[i]]
        matind = z1$ind#RP index and partition index
        eind = which(folds == matind[2])
        enrp[eind, matind[1]] = z1$rpc
        enE[eind, ] = enE[eind, ] + z1$pE1#the sum of the ensemble RP matrix
    }
#     cat("Dim of enrp is", dim(enrp), "\n")
    
    
#     for (j in 1:z[2]) {
# #         cat("The ", j, "-th partition; ", sep = "")
#         # partition
#         eind = which(folds == j)
#         for (i in 1:z[1]) {
# #             cat("The ", i, "-th RP\n", sep = "")
#             # applications of random projection
#             
#             z1 = enlist[[i, j]]
#             if(is.null(z1)){
#                 cat("z1 with", i, "and", j, "is null\n")
#                 cat(str(enlist[i,j]), "\n")
#             }
# #             if(is.null(enlist[i, j])){
# #                 cat("enlist[i, j] with", i, "and", j, "is null\n")
# #             }
# #             cat("The ", j, "-th partition, ", "The ", i, "-th RP: ", names(z1), "\n", sep = "")
#             ptc = z1$rpc
#             pE = z1$pE1
# #             vE = z1$vE1
#             maxsil = z1$maxsil
#                 
# #             z1 = enlist[i, j]
# #             nz = names(z1)
# #             ptc = z1[[nz]]$rpc  #the clustering results by RP for a partitio of the data
# #             pE = z1[[nz]]$pE1  #the partial dim-reduced matrix
# #             maxsil = z1[[nz]]$maxsil
#             
# #             if(sum(is.na(ptc))){#ptc is NULL
# #                 cat("ptc contains NA!\n")
# #             }
#             if(length(eind) == length(ptc)){
#                 enrp[eind, i] = ptc  #the clustering results of different random projection
#             }else{
#                 cat("ptc is", ptc, "and the length of eind is", length(eind), "\n")
#                 enrp[eind, i] = ptc  #the clustering results of different random projection
#             }
#             
#             # print(paste('enrp[folds == j, i] length is: ', length(enrp[folds == j, i]), sep
#             # = ''))
#             
#             # enE[folds == j, ] = enE[folds == j, ] + pE*maxsil
#             
#             enE[eind, ] = enE[eind, ] + pE
#             
#         }
#     }
    
    
    # finalrowColor = wMetaC(enrp) 
    ####parallel programming####
#     frowColor = foreach(t = 1:T, .combine = "c") %dopar% {
    registerDoParallel(n.cores)
    frowColor = foreach(t = 1:T) %dopar% {
        tind = which(folds == t, arr.ind = TRUE)
#         print(paste("the min and max of tind are:", min(tind), " and ", max(tind), sep = ""))
        # tind = seq((t-1)*ng + 1, t*ng) if(t == T){ tind = tind[tind<=ncells] # tind =
        # seq((T-1)*ng + 1, ncells)#select only those within the indices }
        ftmp = wMetaC(enrp[tind, ], hmethod, enpN.cluster, minN.cluster, maxN.cluster, 
            sil.thre, height.Ntimes)  #ensemble of several applications of RPs
        # enmetrics = ARI(gtc[reind[tind], ], ftmp$finalC)#performance evaluation
        # print(enmetrics) ftmp$finalC
        fC = paste(ftmp$finalC, "en", t, sep = "")  #the letters 'en' represent the ensemble clustering
        x0 = ftmp$x0#output matrix for ensemble clustering for further visualization
        
        wres = list()
        wres$fC = fC
        wres$x0 = x0
        wres$nwC = length(unique(fC))#number of clusters
        return(wres)
    }
    stopImplicitCluster()
    
    fColor = character(length = 0)#declare an empty character
    for(t in 1:T){
        fColor = c(fColor, frowColor[[t]]$fC)
    }
	
    uC = unique(fColor)
#     print(uC)
    lenuC = length(uC)#number of meta clusters for sMetaC
    sx0 = matrix(0, ncells, lenuC)
    j1 = 0
    j2 = 0
    for(t in 1:T){
        p = length(frowColor[[t]]$fC)#number of cells in block t
        q = frowColor[[t]]$nwC#number of unique clusters in block t
        sx0[seq(j1+1, j1+p, by = 1), seq(j2+1, j2+q, by = 1)] = frowColor[[t]]$x0##no transpose
#             print(dim(frowColor[[t]]$x0))
#             print(dim(sx0[(j1+1):(j1+p), (j2+1):(j2+q)]))
        j1 = j1 + p
        j2 = j2 + q
    }
    # stop('Trial 1 success!') frowColor = unlist(sapply(1:T, function(t){tind =
    # seq((t-1)*ng + 1, t*ng); if(t == T){tind = seq((T-1)*ng + 1, ncells)}; ftmp =
    # wMetaC(enrp[tind,]); paste(ftmp$finalC, '_', t, sep = '')}))
    
# # #     use fColor instead of frowColor 
    # reorganizing the total cells
    if (T == 1) {
        # if the dataset is not partitioned
#         SrowColor = frowColor
        SrowColor = fColor
        N.cluster = length(unique(SrowColor))  #note that the number of clusters is determined by the unique number in the final round.
        # print(paste("The optimal number of clusters is: ", N.cluster, sep = ""))
        cat("The optimal number of clusters is: ", N.cluster, "\n")
        
        x0 = sx0#for visualization
    } else {
        # newinE = RPmat(data.matrix(E), p)#do the RP;the result is a list E1 =
        # t(newinE$projmat)#for those which need tranpose E1 = data.matrix(E1)
        E1 = enE/K
        
#         stmp = sMetaC(frowColor, E1, folds, hmethod, N.cluster, minN.cluster, 
#             maxN.cluster, sil.thre, height.Ntimes)
        stmp = sMetaC(fColor, E1, folds, hmethod, N.cluster, minN.cluster, 
            maxN.cluster, sil.thre, height.Ntimes, n.cores)
        SrowColor = stmp$finalColor
        # SrowColor = sMetaC(frowColor, E1, folds) SrowColor = sMetaC(frowColor, E, p,
        # folds)
        
        #for visualization
        stf = stmp$tf
#             print(stf)
        sn = length(unique(stf))
        x0 = matrix(0, nrow = nrow(sx0), ncol = sn)
        for(i in 1:sn){
            si = which(stf == i)
            if(length(si)>1){
                x0[, i] = t(rowSums(sx0[, si]))#for visualization
            }else{
                x0[, i] = sx0[, si]
            }
        }
    }
#     finalrowColor = vector(mode = "character", length = length(SrowColor))  #initialization
    finalrowColor = SrowColor
    viE = enE/K
    if(ncol(E) < 1e5){
        finalrowColor[reind] = SrowColor  #reorganizing the final results
    
        x0[reind, ] = x0#reorganizing
    
        viE[reind, ] = viE#reorganizing
    }
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
    
    ####reorganize the final prediction results
    finalrowColor = as.numeric(finalrowColor)
    y = finalrowColor
    uy = unique(y)
    newy = match(y, uy)#we use the index instead
    newuy = unique(newy)
    
    tn = table(newy)
    tn1 = tn[order(as.numeric(names(tn)))]#order the distributions of the clustering results
    
    enresults = list()
    # enresults$finalrowColor = finalrowColor
    enresults$pred_clusters = newy
    enresults$unique_pred_clusters = sort(as.numeric(newuy))
    enresults$distr_pred_clusters = tn1
    # enresults$finalmetrics = finalmetrics
    enresults$N.pred_cluster = length(newuy)
    if(forview){#if we want to visualiza data, we need to save the dimension-reduced matrices/feature vectors
        enresults$x0 = x0
        enresults$viE = viE
    }
    # save(enresults, file=paste(outdir,'enresults_', K, 'times.RData', sep = ''))
    # saveRDS(enresults, file=paste(outdir,'largeresults_', dir1, '.rds', sep = ''))
    return(enresults)
}

#' Test if log-transformation is necessary.
#'
#' Somtimes log-transformation can boost the clustering performance; however, it is not always the case. By this function, SHARP can automatically determine whether a log-transformation is necessary.
#'
#' @param scExp input single-cell expression matrix
#'
#' @param ncells number of cells in the expression matrix
#' 
#' @param p number of reduced dimensions by random projection
#' 
#' @param sncells number of cells randomly selected for checking log-transform. By default, sncells = 100
#' 
#' @examples
#' flag = testlog(scExp)
#'
#' @import Matrix
#'
#' @import foreach
#'
#' @import parallel
#'
#' @import doParallel
#' 
#' @export
testlog<- function(scExp, ncells, p, sncells, n.cores){

    if (missing(sncells)) {
        sncells = 100
    }
    sncells = min(sncells, ncells)#if the dataset is smaller than 100, use the original number of cells
    
    reind = sample(ncells)
    
    E = scExp
    sE = E[, reind[1:sncells]]#selected subset
    #random matrix
    rM0 = ranM(E, p, 5)#reproducible
    rM0 = as.matrix(rM0)
    
    k = 2
    
    registerDoParallel(n.cores)
    msil = foreach(k = 1:2, .combine = c) %dopar% {
    # msil = foreach(k = 1:2, .combine = c) %do% {
        if(k == 2){
            # cat("log-transform when k = ", k, "\n")
            sE = log2(sE + 1)#logarithm transform
        }
        
        newsE = as.matrix(sE)
        sE1 = 1/sqrt(p) * t(rM0) %*% newsE#Random Projection
#         print("works1")
        newsE1 = data.matrix(t(sE1))
    
        tmp = invisible(getrowColor(newsE1, "ward.D", indN.cluster = NULL, 2, 40, 
            0, 2))
        
#         tmp = getrowColor(newsE1, "ward.D", indN.cluster = NULL, 2, 40, 
#                                     0, 2)
        return(tmp$maxsil)
    }
    stopImplicitCluster()
    # msil = tmp$maxsil
#     print(msil)
#     if(msil[1] < msil[2]){#log transform
    if(msil[1] < 0.75 && msil[1] >= 0.95*msil[2]){
        flag = TRUE
    }else{#no need to transform
        flag = FALSE
    }
    return(flag)
}
