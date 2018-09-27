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
SHARP_unlimited2 <- function(scExp, ensize.K, reduced.ndim, partition.ncells, hmethod, N.cluster = NULL, enpN.cluster = NULL, indN.cluster = NULL, minN.cluster, maxN.cluster, sil.thre, height.Ntimes, logflag, n.cores, forview = TRUE, rN.seed) {
    # timing
    start_time <- Sys.time()  #we exclude the time for loading the input matrix
    
    title = "scRNA-Seq Clustering"
    
    if (missing(scExp)) {#scRNA-seq expression matrix
        stop("No expression data is provided!")
    }
    
    ncells = sum(sapply(1:length(scExp), function(i) dim(scExp[[i]])[2]))#number of cells
    ngenes = nrow(scExp[[1]])#number of genes
 
    if (missing(ensize.K)) {
        ensize.K = 5
    }
    
     if (missing(reduced.ndim)) {#reduced dimension
        # default dimensions to be reduced
        reduced.ndim = ceiling(log2(ncells)/(0.2^2))  #reduced 100 times of dimensions; about 200-dim
    }
#     reduced.ndim = ceiling(log2(ncells)/(0.2^2))#the reduced dimension
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
    
    len = length(scExp)#number of partitions/blocks
    
     if (missing(logflag)){#whether we need to check log-transform; 
        if(ncells < 1e4){
            logflag = TRUE
        }else{#if the number of cells is too large, we do not check but directly adopt the log-transform
            logflag = FALSE
        }
        
    }
    
    if (logflag){#if we want to check whether log-transform or not
        #by default, we randomly select 100 cells to determine log transform or not
        sncells = 100
   
        nc1 = ncol(scExp[[1]])
        flag = testlog(scExp[[1]], nc1, reduced.ndim, sncells)#test whether log transform is necessary or not; by default, it will randomly select 100 cells to test
        
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
    
     if (missing(n.cores)) {
        # number of cores to be used, the default is to use all but one cores
        n.cores = detectCores() - 1
    }
    # registerDoMC(n.cores)
    registerDoParallel(n.cores)
    
  
    
    
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
    
    
    p = reduced.ndim         
    rM = foreach(k = 1:ensize.K) %dopar% {###to use consistent random matrices across different partitions
            if (rN.seed == 0.5) {
                ranM2(ngenes, p, 0.5)  #for random results
            } else {
                ranM2(ngenes, p, 50 + rN.seed + k)  #for reproducible result
            }
    }
    
   
        
    ifColor = vector("list", length = len)
    iE1 = vector("list", length = len)
    ifolds = vector("list", length = len)
    y = list()
    k = 0 
    for(i in 1:len){#for each partition file
        cat("=======================================================================\nProcessing Partition", i, "of the scRNA-seq data...\n")
#         mat = readRDS(allfiles[i])
        y = SHARP_fpart(scExp[[i]], ensize.K, reduced.ndim, partition.ncells, hmethod, N.cluster, enpN.cluster, indN.cluster, minN.cluster, maxN.cluster, sil.thre, height.Ntimes, flag, rM, rN.seed)
    
        
#         mat = as.matrix(scExp[[i]])
#         cat("=======================================================================\n")
        
        k = k + y$ncells
        ifolds[[i]] = y$folds + k
        ifColor[[i]] = paste(y$fColor, "s", i, sep = "")
        iE1[[i]] = t(y$E1)###for ease of ordering; by default, R concatenate elements by columns; but here the feature dim (column) is the same for each element of the list
        ngenes = y$ngenes
        
        rm(y)
        gc()
    }
    
#     ncells = k
    
    
    folds = unlist(ifolds)#convert a list to a vector
    fColor = unlist(ifColor)
    E1 = t(matrix(unlist(iE1), ncol = ncells))
    rm(iE1)
    gc()
    
    cat("=======================================================================\n")
    cat("Length of fColor:", length(fColor), "\n")
    cat("Number of unique elements of fColor:", length(unique(fColor)), "\n")
    cat("Dim of E1:", dim(E1), "\n")
#     print(y[[1]]$paras)
    
    cat("=======================================================================\n")
    cat("Meta-clustering for combining blockwise results...\n")
    stmp = sMetaC(fColor, E1, folds, hmethod, N.cluster, minN.cluster, 
            maxN.cluster, sil.thre, height.Ntimes)
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
    
    enresults = list()
    enresults$pred_clusters = finalrowColor
    enresults$unique_pred_clusters = sort(as.numeric(uf))
    enresults$distr_pred_clusters = tn1
    enresults$N.pred_clusters = luf
    enresults$N.cells = ncells
    enresults$N.genes = ngenes
    enresults$reduced.ndim = reduced.ndim
    enresults$ensize.K = ensize.K
    if(forview){
        enresults$viE = E1
        
#         x0 = matrix(0, ncells, luf)
#         for(i in 1:ncells){
#             x0[i, finalrowColor[i] == uf] = 1
#         }
        x0 = sparseMatrix(i = 1:ncells, j = as.numeric(finalrowColor), x = 1, dims = c(ncells, luf))
        enresults$x0 = x0
    }
    enresults$time = t
    
     ###parameter saving###
    paras = list()
    paras$ensize.K = ensize.K
    paras$reduced.ndim = reduced.ndim
    paras$partition.ncells = partition.ncells
    paras$hmethod = hmethod
    paras$N.cluster = N.cluster
    paras$minN.cluster = minN.cluster
    paras$maxN.cluster = maxN.cluster
    paras$sil.thre = sil.thre
    paras$height.Ntimes = height.Ntimes
    paras$n.cores = n.cores
    
    enresults$paras = paras
    
    stopImplicitCluster()
    
    
    
    return(enresults)
}

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
SHARP_fpart <- function(scExp, ensize.K, reduced.ndim, partition.ncells, hmethod, 
    N.cluster, enpN.cluster, indN.cluster, minN.cluster, maxN.cluster, sil.thre, 
    height.Ntimes, flag, rM, rN.seed){

    start_time <- Sys.time()  #we exclude the time for loading the input matrix
    
    ncells = ncol(scExp)#number of cells
    ngenes = nrow(scExp)
    cat("Number of cells:", ncells, "\n")
    cat("Number of genes:", ngenes, "\n")
     ######## Partition the large data into several groups###########
    p = reduced.ndim
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
        rM = foreach(k = 1:K) %dopar% {
            if (rN.seed == 0.5) {
                ranM(E, p, 0.5)  #for random results
            } else {
                ranM(E, p, 50 + rN.seed + k)  #for reproducible result
            }
        }
    }
    
    
    #### parallel programming#### rerowColor = foreach(t=1:T, .combine = 'c') %dopar%{
    enlist = foreach(k = 1:K, .combine= "c", .packages= "Matrix") %:% foreach(t = 1:T) %dopar% {
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
#         E1 = 1/sqrt(p) * t(rM[[k]]) %*% newE
        # newinE = RPmat(inE, p)#do the RP;the result is a list E1 =
        # t(newinE$projmat)#for those which need tranpose
        
        E1 = t(E1)
#         newE1 = t(E1)
        
        newE1 = data.matrix(E1)  #convert the spase dim-reduced matrix back to full dim-reduced matrix
        newE1 = round(newE1, digits = 1)

#         rm(list = c(newE, inE, E1))
#         gc()

        # if(ncells>10000){#for large-scale datasets, we normalize the values for
        # memory-efficient computation E1 = round(E1/100, digits = 3)#due to the high
        # dimensions (>10000) in the original data, the values of E1 may be very large,
        # better to normalize to a small scale for later clustering }
#         cat("Process ", pid, "----", Sys.time(), "\n", sep = "")
#         cat("Dim reduction completed for:", k, "-th RP and", t, "-th partition\n")
        maxN.cluster = 40#for partition clustering
        tmp = getrowColor(newE1, hmethod, indN.cluster, minN.cluster, maxN.cluster, 
            sil.thre, height.Ntimes)  #hierarchical clustering; the result is a list containing both the predicted clusters and the maximum silhouette index
#         cat("Process ", pid, "----", Sys.time(), "\n", sep = "")
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
        
#         rm(newE1)
#         rm(list = c(newE, inE, E1, newE1))
        gc()
        
#         cat(str(tmplist), "for:", k, "-th RP and", t, "-th partition\n")
        
#         if(length(tmplist) == 0){
#             cat("tmplist is null!\n")
#         }
        return(tmplist)
        # rerowColor[rind] = paste(tmp, '_', t, sep = '')#for distinguishing for each
        # smaller group clustering
    }
    
    
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
#         cat("enE[eind, ]:", dim(enE[eind, ]), "\n")
#         cat("z1$pE1", dim(z1$pE1), "\n")
        enE[eind, ] = enE[eind, ] + z1$pE1#the sum of the ensemble RP matrix
    }
#     cat("Dim of enrp is", dim(enrp), "\n")
    
      # finalrowColor = wMetaC(enrp) 
    ####parallel programming####
#     frowColor = foreach(t = 1:T, .combine = "c") %dopar% {
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
    
    
    E1 = enE/K
    
    ###readjust back
    if(ncol(E) < 1e5){
        fColor[reind] = fColor
        E1[reind, ] = E1
        folds[reind] = folds
    }
    
    nmcluster = length(unique(fColor))
    cat("Number of meta clusters:", nmcluster, "\n")
    
    #save results
    res = list()
    res$fColor = fColor
    res$E1 = E1
    res$nmcluster = nmcluster
    res$folds = folds
    res$ncells = ncells
    res$ngenes = ngenes
    
    end_time <- Sys.time()
    
    t <- difftime(end_time, start_time, units = "mins")  #difference time in minutes
    cat("Running time:", t, "minutes\n")
    return(res)
}
