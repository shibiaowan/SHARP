#' similarity-based ensemble meta-clustering
#'
#' This function is to do similarity-based ensemble clustering for combining the results from different partitions of single cells. The similarity is measured by the group-group correlations.
#'
#' @param nC a m*n matrix, where m is the number of cells, n is the number of clustering algorithms, and the (i,j)-element of nC represents the cluter for the i-th cell by the j-th clutering predictor.
#'
#' @examples
#' finalrowColor = sMetaC(rerowColor, sE1, folds)
#'
#' @import cluster
#'
#' @import clues
#'
#' @import clusterCrit
#'
#' @export
sMetaC <- function(rerowColor, sE1, folds, hmethod, finalN.cluster, minN.cluster, 
    maxN.cluster, sil.thre, height.Ntimes) {
    # This is to do meta-clustering the results obtained for each smaller group of
    # the original large-scale datasets
    R = unique(rerowColor)  #unique clusters
    # print(R)
    # cat("Unique clusters for similarity-based meta-clustering are: \n", R, "\n")
    
    nC = length(unique(rerowColor))  #number of unique clusters
    
    # get the number of clusters for each partition groups
    T = unique(folds)
    t0 = vector(mode = "numeric", length = length(T))
    # ff = numeric()
    for (t in 1:length(T)) {
        tmp = which(folds == T[t])
        t0[t] = length(unique(rerowColor[tmp]))  #number of clusters for each partition
        # ff = c(ff, rep(t, each = t0[t]))
    }
    
    
    # pind = which.max(t0)#the partition which has the most predicted clusters rind =
    # setdiff(1:T, pind) #preprocessing
    sE1 = sE1[apply(sE1, 1, function(x) sd(x) != 0), ]  #avoid NA when sd is 0 during scaling
    sE1 <- t(scale(t(sE1)))  #vector scaling
    
    # mE = matrix(0, nrow = length(R), ncol = dim(sE1)[2]) for(i in 1:length(R)){ G =
    # sE1[rerowColor == R[i], , drop = FALSE] mE[i, ] = colMeans(G)*dim(G)[1] } mf =
    # getwfperC(frowColor, E, p)#weighted features for each partitioned cluster
    cb = combn(nC, 2)  #all possible combinations (n*(n-1)/2)
    # t = p alls = apply(cb, 2, getpaircor, rerowColor = rerowColor, E1 =
    # E1)#calculate the correlation/simialrity for all combinations (apply to the
    # columns) parallel programming####
    alls = foreach(t = 1:dim(cb)[2], .combine = "c") %dopar% {
        getpaircor(cb[, t], rerowColor, sE1, R)
    }
    # alls = foreach(t=1:dim(cb)[2], .combine = 'c') %dopar%{exp(-sum((mf[cb[1, t]] -
    # mf[cb[2, t]])^2)/t)}
    
    S0 = sparseMatrix(i = cb[1, ], j = cb[2, ], x = alls, dims = c(nC, nC))  #triangle part of the S
    S = S0 + t(S0) + diag(nC)
    
    # w = matrix(0, nrow = length(R), ncol = 2) for(i in 1:length(R)){ tind
    # =setdiff(1:length(R), which(ff == ff[i]))#rest index w[i,1] = min(S[i, tind]) }
    
    # # S = t(t(S)/(colSums(S)-1)) # diag(S) = 1 print(S)
    
    
    # S = matrix(0, ncol = nC, nrow = nC)#initialization for a similarity matrix of
    # all the clusters for(i in 1:nC){ for(j in (i+1):nC){ G1 = E1[rerowColor ==
    # allC[i], ]#find all the cells for cluster i G2 = E1[rerowColor == allC[j],
    # ]#find all the cells for cluster j cor1 = cor(t(G1), t(G2))#correlation between
    # these two clusters S[i, j] = (max.col(cor1) + max.col(t(cor1)))/(nrow(G1) +
    # nrow(G2))#similarity between two clusters S[j, i] = S[i, j] } S[i, i] = 1 }
    
    # res = kmeans(mE, max(t0)) d=as.dist(1-cor(t(mE)))
    if (missing(minN.cluster)) {
        minN.cluster = max(2, min(t0))
    }
    hres = get_opt_hclust(S, hmethod, N.cluster = finalN.cluster, minN.cluster, maxN.cluster, 
        sil.thre, height.Ntimes)
    
    tf = hres$f
    v = hres$v
    # d=as.dist(1-S) # d=as.dist(S) # print(S) # metah = hclust(d, method='ward.D') h
    # = hclust(d, method='ward.D')#although it is also meta, it is for a single large
    # dataset # S = mE # nc = 2:min(40, dim(S)[1]-1) nc = max(2, min(t0)):min(40,
    # dim(S)[1]-1) v = matrix(0, nrow = dim(S)[1], ncol = length(nc))#for all
    # different numbers of clusters msil = rep(0, length(nc))#declare a vector of
    # zeros CHind = rep(0, length(nc)) print(paste('The height for the top 10 are: ',
    # tail(h$height, n = 10), sep = '')) # stop('Stop here!') for(i in 1:length(nc)){
    # # foreach(i=1:length(nc)) %dopar%{ v[,i] = cutree(h, k = nc[i])#for different
    # numbers of clusters sil = silhouette(v[,i], d)#calculate the silhouette index #
    # msil[i] = mean(sil[,3])#the mean value of the index msil[i] =
    # median(sil[,3])#the mean value of the index # CHind[i] = get_CH(S, v[,i],
    # disMethod = 'Euclidean') ch0 =
    # intCriteria(data.matrix(S),as.integer(v[,i]),'Calinski_Harabasz') CHind[i] =
    # unlist(ch0, use.names = F)#convert a list to a vector/value } print(msil)#the
    # average sil index for all cases print(CHind) # oind = which.max(msil)#the index
    # corresponding to the max sil index tmp = which(msil == max(msil))#in case there
    # are more than one maximum if(length(tmp)>1){oind =
    # tmp[ceiling(length(tmp)/2)]}else{oind = tmp} if(max(msil)<=0.35){ oind =
    # which.max(CHind) if(oind ==1){#it's likely that the CH index is not reliable
    # either tmp = tail(h$height, n =10)#the height diftmp = diff(tmp) flag = diftmp
    # > tmp[1:(length(tmp)-1)]#require the height is more than 2 times of the
    # immediate consecutive one if(any(flag)){#if any satifies the condition; make
    # sure at least one satisfied pind = which.max(flag) opth = (tmp[pind] +
    # tmp[pind+1])/2#the optimal height to cut optv = cutree(h, h = opth)#using the
    # appropriate height to cut oind = length(unique(optv)) - 1#for consistency } } }
    # # oind = 1 tf = v[,oind]#the optimal clustering results
    
    # tf=cutree(metah,k=N.cluster)
    rerowColor[] <- vapply(rerowColor, function(x) tf[match(x, R)], numeric(1))  #apply to every element;this is for matrices
    
    # finalC = apply(rerowColor, 1, function(d)
    # names(sort(table(d),decreasing=TRUE)[1])) #find the most repeated elements for
    # each row
    
    # f = rerowColor my = sE1 silthre = 0 nthre = 0.05*length(f) nf =
    # character(length = length(f)) # if(hres$maxsil > silthre){ # nf =
    # as.character(f) # }else{ uf = unique(f) r = length(unique(f)) for(t in 1:r){
    # tind = which(f==uf[t]) if(length(tind)<= nthre){ nf[tind] =
    # as.character(f[tind]) }else{ tE = my[tind,] d1=as.dist(1-cor(t(tE))) hres1 =
    # gethclust(d1, tE) rf = hres1$f nf[tind] = paste(rf, 'sub', t, sep = '') } } # }
    
    
    finalColor = rerowColor
    
    N.cluster = length(unique(finalColor))  #note that the number of clusters for meta-clustering is not determined by previous selection, but by the unique number in the final round.
    
    # print(paste('The optimal number of clusters for combining partitions is: ',
    # N.cluster, sep = ''))
    cat("The optimal number of clusters for combining partitions is:", N.cluster, 
        "\n")
    
    finalres = list()
    finalres$finalColor = finalColor
    finalres$tf = tf#optimal clustering results
    return(finalres)
    # return(finalColor)
}



getpaircor <- function(pairind, rerowColor, sE1, R) {
    G1 = sE1[rerowColor == R[pairind[1]], , drop = FALSE]  #find all the cells for cluster i; 'drop = F' is to avoid a one-row/column matrix being converted to a vector
    G2 = sE1[rerowColor == R[pairind[2]], , drop = FALSE]  #find all the cells for cluster j
    
    # cor1 = cor(t(G1), t(G2))#correlation between these two clusters t =
    # 1.5*dim(G1)[2]
    aG1 = colMeans(G1)
    aG2 = colMeans(G2)
    # corss = dist(rbind(aG1, aG2)) corss = exp(-sum((aG1 - aG2)^2)/t)
    corss = cor(aG1, aG2)
    # as1 = 1/(1 + exp(-aG1)) as2 = 1/(1 + exp(-aG2)) ass =
    # sum(aG1*aG2)/sqrt(sum(aG1^2)*sum(aG2^2)) corss = 1/(1 + exp(-ass)) cor1 =
    # foreach(i = 1:dim(G1)[1], .combine = 'cbind') %:% foreach(j = 1:dim(G2)[1],
    # .combine = 'c') %do%{ # sum(G1[i,]*G2[j,])/sqrt(sum(G1[i,]^2)*sum(G2[j,]^2))
    # exp(-sum((G1[i,] - G2[j,])^2)/t) } m = dim(cor1) if(m[1] <= m[2]){#use the
    # smaller-size group to measure the correlation/similarity corss =
    # sum(do.call(pmax, data.frame(cor1)))/nrow(G1)#similarity between two clusters;
    # summing all the max correlation/similarity }else{ corss = sum(do.call(pmax,
    # data.frame(t(cor1))))/nrow(G2) }
    
    # c1 = as.vector(cor1) c2 = c1[order(-c1)]#sorting the number in descending order
    # corss = sum(c2[1:ceiling(length(c2)/2)])/sum(c2) c1 = do.call(pmax,
    # data.frame(cor1)) c2 = do.call(pmax, data.frame(t(cor1))) i0 = intersect(c1,
    # c2) u0 = union(c1, c2) corss = sum(i0)/sum(u0) corss = (sum(do.call(pmax,
    # data.frame(cor1))) + sum(do.call(pmax, data.frame(t(cor1)))))/(nrow(G1) +
    # nrow(G2))#similarity between two clusters; summing all the max
    # correlation/similarity corss = (sum(do.call(pmax, data.frame(cor1)))/nrow(G1) +
    # sum(do.call(pmax, data.frame(t(cor1))))/nrow(G2))/2 corss = max(max(cor1))
    # print(corss) corss = mean(mean(cor1))
    return(corss)
}
