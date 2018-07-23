#' weighted ensemble clustering
#'
#' This function is to do weighted ensemble clustering for meta-clustering.
#'
#' @param nC a m*n matrix, where m is the number of cells, n is the number of clustering algorithms, and the (i,j)-element of nC represents the cluter for the i-th cell by the j-th clutering predictor.
#'
#' @examples
#' finalrowColor = wMetaC(enrp)
#'
#' @import cluster
#'
#' @import clues
#'
#' @export
wMetaC <- function(nC, hmethod, enN.cluster, minN.cluster, maxN.cluster, sil.thre, 
    height.Ntimes) {
    # This is to obtain the weight matrix for each cluster solution for following
    # meta-clustering
    N = nrow(nC)  #number of points
    C = ncol(nC)  #number of clustering methods/times; or K
    # AA = Matrix(0, nrow = N, ncol = N, sparse = T)#initialization for (i in 1:C){
    # AA = AA + getA(nC[,i]) }
    
    AA = Reduce("+", apply(nC, 2, getA))  #sum over obtained matrices; execute along the column and then matrix sum
    AA = AA/C
    
    # vtmp = prcomp(AA, rank. = 2)#select the top-2 PCs vAA = vtmp$x#the 2-dim data
    
    # AA = sum(apply(nC, 2, getA))/C
    indA = which(AA != 0, arr.ind = T)  #find non-zero indices of AA
    nd = vapply(AA[indA], function(x) x * (1 - x), numeric(1))
    # ####parallel programming#### nd = foreach(t=1:dim(indA)[1], .combine = 'c')
    # %dopar%{AA[indA[t,]]*(1-AA[indA[t,]])}#calculate the weight s for all
    # combinations
    
    newAA = sparseMatrix(i = indA[, 1], j = indA[, 2], x = nd, dims = c(N, N))
    
    
    
    # AA[] <- vapply(AA, function(x) x*(1-x), numeric(1))#apply to every element
    w0 = 4/N * rowSums(newAA)  #the weight for each point
    # e = 0.001
    e = 0.01
    w1 = (w0 + e)/(1 + e)  #adjusted point weight
    
    # #revise the labels for differentiation x <- vector(mode='character', length=0)
    # for (i in 1:C){ # t1 = factor(nC[,i])#rowColor is a vector like 'red red purple
    # brown ...'  # levels(t1) = paste(levels(t1), '_', i, sep = '')#make each
    # cluster in each column different t1 = paste(nC[,i], '_', i, sep = '') x = c(x,
    # t1)#reshape the matrix to a vector } #Cluster weight value newnC <- matrix(x,
    # nrow = N, byrow = FALSE)#reshape a vector to a matrix; by column
    
    # revise the labels for differentiation # newnC = matrix('0', nrow = N, ncol =
    # C)#revised cluster name x = character(length = 0)#initialization for (i in
    # 1:C){ newnC = paste(nC[,i], '_', i, sep = '')#distinguish the cluster labels
    # for each column/solution x = c(x, newnC)#extend the matrix to a vector for ease
    # of processing # newnC[,i] = paste(nC[,i], '_', i, sep = '')#distinguish the
    # cluster labels for each column/solution # x = c(x, newnC[,i])#extend the matrix
    # to a vector for ease of processing }
    x = as.vector(sapply(1:C, function(i) {
        paste(nC[, i], "_", i, sep = "")
    }))  #convert the matrix (N*C) to vector (concatenating them)
    
    newnC <- matrix(x, nrow = N, byrow = FALSE)  #reshape a vector to a matrix; by column
    
    R = unique(x)  #all unique labels
    allC = length(R)  #number of all unique labels
    # S = matrix(0, allC, allC)
    
    cb = combn(allC, 2)  #all possible combinations (n*(n-1)/2)
    alls = apply(cb, 2, getss, R = R, x = x, w1 = w1)  #calculate the weight s for all combinations
    
    # ####parallel programming#### alls = foreach(t=1:dim(cb)[2], .combine = 'c')
    # %dopar%{getss(cb[,t], R, x, w1)}#calculate the weight s for all combinations
    
    S0 = sparseMatrix(i = cb[1, ], j = cb[2, ], x = alls, dims = c(allC, allC))  #triangle part of the S
    S = S0 + t(S0) + diag(allC)
    # for (k in c(1:(allC-1))){ k1 = which(x %in% R[k])#find samples with k-th
    # cluster d1 = as.numeric(unlist(strsplit(R[k], '_'))[2])#the name contains only
    # two parts; get the numbering part #convert d1 to an integer newk1 = k1 - (d1 -
    # 1)*N#the index for (j in c((k + 1):allC)){ k2 = which(x %in% R[j]) #%%
    # N#modulus over N d2 = as.numeric(unlist(strsplit(R[j], '_'))[2])#the name
    # contains only two parts; get the numbering part newk2 = k2 - (d2 - 1)*N#the
    # index intset = intersect(newk1, newk2)#set intersection uset =
    # union(newk1,newk2)#set union if (length(intset) != 0){#if 0, no need to
    # calculate S[k, j] = sum(w1[intset])/sum(w1[uset]) S[j, k] = S[k, j]#symmetric }
    # } S[k, k] = 1 }
    
    
    # if ((length(k1) != 0) && (length(k2) != 0)){ # if (k1 == 0){ # k1 = N # } k1[k1
    # == 0] = N#the final one k2[k2 == 0] = N # if (k2 == 0){ # k2 = N # } intset =
    # intersect(k1, k2) uset = union(k1,k2) if (length(intset) != 0){ S[k, j] =
    # sum(w1[intset])/sum(w1[uset]) # S[j, k] = S[k, j]#symmetric } } }
    if (missing(sil.thre)) {
        # if sil.thre is not assigned
        sil.thre = 0
    }
    hres = get_opt_hclust(S, hmethod, N.cluster = enN.cluster, minN.cluster, maxN.cluster, 
        sil.thre, height.Ntimes)  #solely using the silhouette index as the criteria
    
    tf = hres$f
    v = hres$v
    cat("The number of clusters before voting is: ", hres$optN.cluster, "\n")
    
    # d=as.dist(1-S) # print(S) # metah = hclust(d, method='ward.D') metah =
    # hclust(d, method='ward.D') # print(dim(S)) maxc = min(40, dim(S)[1]-1)#maximum
    # number of clusters; at most (maxc-1) clusters because maxc clusters mean no
    # needing to cluster v = matrix(0, nrow = dim(S)[1], ncol = maxc-1)#for all
    # different numbers of clusters # print(dim(v)) msil = rep(0, maxc-1)#declare a
    # vector of zeros # mdunn = rep(0, maxc-1)#for dunn index # mdb = rep(0,
    # maxc-1)#for DB index # wss = rep(0, maxc-2)#within-cluster sum of squares #
    # CHind = rep(0, maxc-1)#CH index allc = c(2:maxc) for(i in 1:length(allc)){ v[,
    # i] = cutree(metah, k = allc[i])#for different numbers of clusters #
    # print(v[,i]) sil = silhouette(v[, i], d)#calculate the silhouette index #
    # msil[i-1] = mean(sil[,3])#the mean value of the index msil[i] =
    # median(sil[,3])#the mean value of the index # mdunn[i-1] = dunn(d, v[,i-1]) #
    # db = index.DB(d, cl = v[, i-1]) # mdb[i-1] = db$DB # spl <- split(d, v[, i-2])
    # # wss[i-2] <- sum(sapply(spl, wss))#within-cluster sum of squares # print('OK')
    # # CHind[i] = get_CH(S, v[, i], disMethod = '1-corr') # print(CHind[i]) } # oind
    # = which.max(msil)#the index corresponding to the max sil index tmp = which(msil
    # == max(msil))#in case there are more than one maximum if(length(tmp)>1){oind =
    # tmp[ceiling(length(tmp)/2)]}else{oind = tmp} print(msil)#the average sil index
    # for different numbers of clusters # print('OK!') # print(CHind) #
    # if(max(msil)<=0.2){ # difCH = diff(CHind) # x1 = which(difCH<0) #
    # if(length(x1)>0){ # oind = min(x1)#local maximum # }else{ # oind = 1 # } # #
    # oind = which.max(CH) # } # #if the silhouette index is not reliable, we use the
    # gap statistics # if(max(msil)<=0.1){ # # # gskmn <- clusGap(S, FUN = kmeans,
    # nstart = 20, K.max = min(maxc, 15), B = min(50, nrow(S))) # g = gskmn$Tab # gap
    # = g[, 'gap']#the gap info # print(gap) # # # oind = which.max(gap)#maximum gap
    # # # sesim = g[, 'SE.sim']#standard error of the gap # print(sesim) # oind =
    # maxSE(gap, sesim)#maximize the gap with parsimony of the model # # # if(oind
    # >=floor(maxc*0.8)){#if the gap stastic keeps increasing until very big number
    # of cluster, we use the first SE-rule (Gap(k) >= Gap(k+1) - SE) # # sesim = g[,
    # 'SE.sim']#standard error of the gap # # print(sesim) # # oind = maxSE(gap,
    # sesim)#maximize the gap with parsimony of the model # # } # } # print(wss) #
    # print(mdunn) # print(mdb) # oind = 3 tf = v[,oind]#the optimal clustering
    # results
    
    
    # tf=cutree(metah,k=N.cluster)
    newnC[] <- vapply(newnC, function(x) tf[match(x, R)], numeric(1))  #apply to every element; reorganizing the clusters for different results
    
    finalC = apply(newnC, 1, function(d) names(sort(table(d), decreasing = TRUE)[1]))  #find the most repeated elements for each row
    
    N.cluster = length(unique(finalC))  #note that the number of clusters for meta-clustering is not determined by previous selection, but by the unique number in the final round.
    
    # print(paste('The optimal number of clusters for ensemble clustering is: ',
    # N.cluster, sep = ''))
    cat("The optimal number of clusters for ensemble clustering is:", N.cluster, 
        "\n")
    
    
    # #select representative cells/data ncp = 100#number of cells per cluster
    # selected seind = numeric()#empty/initialization s_newnC = apply(newnC, 1,
    # sort(table(d),decreasing=TRUE)[1]/length(d))#the score of the most-voting one
    # for(j in 1:N.cluster){ nj = which(finalC == j)#number of cells/data in cluster
    # j tmp = s_newnC[nj]#selected scores njo = order(-tmp)#order in descending order
    # if(length(njo) >= ncp){#if larger than threhold, we selected the first 100
    # data/cells njj = nj[njo[1:ncp]] seind = c(seind, njj) }else{#else just select
    # all seind = c(seind, nj) } } f_seind = finalC[seind]#the predicted clusters for
    # the selected cells
    
    
    # For ease of visualization
    uC = unique(finalC)#unique clusters
# print(uC)
# #   x0 = apply(newnC, 1, function(x){t = rep(0, N.cluster); for(i in c(1:N.cluster)){t[i] = length(which(x %in% i))}; return(t)})
  y0 = apply(newnC, 1, function(x){t = rep(0, N.cluster); for(i in c(1:N.cluster)){t[i] = length(which(x %in% uC[i]))}; return(t)})#need to reorganize before counting
#   print(dim(y0))
  y0 = t(y0)#transpose
  
  
  x0 = matrix(0, nrow = N, ncol = N.cluster)
#   print(dim(x0))
  
  
  tw = 0.5
#   print(uC)
  for(i in 1:N){
    xind = which(finalC[i]==uC)
    x0[i, xind] = 1#the correct clustering result
    allind = which(y0[i,]!=0)#all the counts
    diffind = setdiff(allind, xind)#some other counts which are not the correct cluster
    if(length(diffind) != 0){
        x0[i, diffind] = tw* y0[i, diffind]/y0[i, xind]#use a reduced weight
    }
  }
  
  
#     x0 = apply(newnC, 1, function(x) {
#         t = rep(0, N.cluster)
#         for (i in c(1:N.cluster)) {
#             t[i] = length(which(x %in% i))
#         }
#         return(t)
#     })
#     x0 = t(x0)  #transpose
    
    out = list()  #declare
    out$finalC = finalC
    out$x0 = x0
    # out$seind = seind out$f_seind = f_seind out$pcaAA = vAA out$AA = AA out$S = S
    # return(finalC)
    return(out)
}



#' Obtain the co-location matrix
#'
#' This function is to obtain the weighted co-location matrix from the clustering result for following ensemble clustering.
#'
#' @param rowColor the clustering results by some clustering algorithms
#'
#' @examples
#' AA = getA(rowColor)
#'
#' @import Matrix
#'
#' @export
getA <- function(rowColor) {
    # This is to obtain the weighted co-association matrix for clustering solution
    # rowColor
    N = length(rowColor)  #number of points
    
    # ------------------------------------------------------- #convert a
    # (categorical) vector to a numeric vector t1 = factor(rowColor)#rowColor is a
    # vector like 'red red purple brown ...'  levels(t1) =
    # c(1:length(levels(t1)))#convert the categorical to numeric t1 =
    # as.numeric(as.character(t1))
    # -------------------------------------------------------
    
    # t1 = factor(rowColor)#rowColor is a vector like 'red red purple brown ...'
    # levels(t1) = c(1:length(levels(t1)))#convert the categorical to numeric t1 =
    # as.numeric(as.character(t1))
    
    L = levels(factor(rowColor))
    # A = matrix(0, N, N)#declare a 0-value matrix A = Matrix(0, nrow = N, ncol = N,
    # sparse = T)#use a parse matrix for efficient storage ind = which(rowColor[, c]
    # %in% k)#the same cluster
    
    # find indices for each cluster, then all combinations of indices
    tmp = sapply(L, function(k) {
        r = which(rowColor %in% k)
        expand.grid(r, r)
    })
    # ####parallel programming#### tmp = foreach(t=1:length(L), .combine = 'c')
    # %dopar%{}#calculate the weight s for all combinations
    
    # reshape to the indices
    allind = matrix(unlist(t(tmp)), ncol = 2, byrow = F)  #need transpose
    A = sparseMatrix(i = allind[, 1], j = allind[, 2], x = 1, dims = c(N, N))  #non-zero entries
    
    # for (k in L) { r = which(rowColor %in% k)#for the data in the same cluster with
    # value k ind = expand.grid(r,r)#all the potential combinations tmp =
    # sparseMatrix(i = ind[,1], j = ind[,2], x= 1, dims = c(N,N))#non-zero entries A
    # = A + tmp # A[which(rowColor %in% k), which(rowColor %in% k)] = 1 }
    
    # r = which(rowColor %in% k); ind = expand.grid(r,r) apply(L, function(k)
    # A[which(rowColor %in% k), which(rowColor %in% k)] = 1)
    return(A)
}



#' Obtain the co-location matrix
#'
#' This function is to obtain the weighted co-location matrix from the clustering result for following ensemble clustering.
#'
#' @param rowColor the clustering results by some clustering algorithms
#'
#' @examples
#'  x = as.vector(sapply(1:C, function(i){paste(nC[,i], '_', i, sep = '')}))#convert the matrix (N*C) to vector (concatenating them)
#'  R = unique(x)#all unique labels
#'  alls = apply(cb, 2, getss, R = R, x = x, w1 = w1)#calculate the weight s for all combinations
#'
#' @export
getss <- function(pind, R, x, w1) {
    # This is to get the element of S
    pairk = lapply(pind, getnewk, R = R, x = x, N = length(w1))  #run for two indices
    
    intset = intersect(unlist(pairk[1]), unlist(pairk[2]))  #set intersection
    
    ss = 0
    if (length(intset) != 0) {
        uset = union(unlist(pairk[1]), unlist(pairk[2]))  #set union
        ss = sum(w1[intset])/sum(w1[uset])
    }
    return(ss)
}

getnewk <- function(k, R, x, N) {
    # This is to get the original index of the sample
    k1 = which(x %in% R[k])  #find samples with k-th cluster
    d1 = unlist(strsplit(R[k], "_"))  #the name contains only two parts; get the numbering part
    d = as.numeric(tail(d1, n = 1))  #the last element of the split arrays
    newk1 = k1 - (d - 1) * N  #the index
    return(newk1)
}
