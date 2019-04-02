#' get the optimal hierarchical clustering results with the optimal number of clusters
#'
#' This function is to estimate the optimal number of clusters by combining three indices, including Silhouette index, Calinski-Harabasz (CH) index and height difference.
#'
#' @param mat   either a feature matrix or a similarity matrix derived from the single-cell expression matrix
#'
#' @param hmethod   agglomeration method for hierarchical clustering, the default is 'ward.D'. Certainly, some other methods can also be used, like 'ward.D2', 'single', 'complete', 'average' (= UPGMA), 'mcquitty' (= WPGMA), 'median' (= WPGMC) or 'centroid' (= UPGMC).
#'
#' @param minN.cluster  minimum number of clusters to be tested, the default is 2.
# 
#' @param maxN.cluster  maximum number of clusters to be tested, the default is 40 or equal to the number of cells (within the specific clustering problem) minus 1, whichever is smaller.
# 
#' @param sil.thre  the threshold for the maximum Silhouette index (msil), the default is 0.35. If msil < sil.thre, we should use CH index.
# 
#' @param height.Ntimes  the threshold for the height difference between two adjacent descending-ordered heights obtained after hierarchical clustering. If the height difference is above the threshold, we cut at the median height between the first height and its immediate next which satisfy the criteria.
# 
#' @details Specifically, we first select the maximum Silhouette index (msil) as the reference. If msil > threshold (here we use sil.thre as the threshold, the default value is 0.35), then its corresponding number of clusters is the optimal; otherwise, we use the maximum CH index as the reference. If the number of clusters with the maximum CH index is not 2, then it is the optimal number of clusters; otherwise, we use the adjacent height difference (which is derived from hierarchical clustering). If the former height is larger than a threshold (Ntimes larger than the immediate latter height), then we cut at the mean height between these two, and the corresponding number of clusters is the optimal one; otherwise, we do not cut.
#'
#' @return a list containing the optimal hierarchical clustering results, the optimal number of clusters, the corresponding maximum Silhouette index and other indices.
#'
#' @examples
#' hres = get_opt_hclust(mat)
#'
#' @import cluster
#'
#' @import clues
#'
#' @import clusterCrit
#'
#' @export
get_opt_hclust <- function(mat, hmethod, N.cluster, minN.cluster, maxN.cluster, sil.thre, 
    height.Ntimes) {
    # if no agglomeration method for hierarchical clustering is provided
    if (missing(hmethod) || is.null(hmethod)) {
        hmethod = "ward.D"  #the default hierarchical clustering agglomeration method is 'ward.D'
    }
    
    # if no minimum number of clusters is provided
    if (missing(minN.cluster) || is.null(minN.cluster)) {
        minN.cluster = 2  #by default, we try the minimum number of clusters starting from 2
    }
    
    # if no maximum number of clusters is provided
    if (missing(maxN.cluster) || is.null(maxN.cluster)) {
        maxN.cluster = 40  #by default, we try the maximum number of clusters as large as 40 or the number of cells minus 1, whichever is smaller.
    }
    
    # if no threshold for the maximum Silhouette index is provided
    if (missing(sil.thre) || is.null(sil.thre)) {
        sil.thre = 0.35  #by default, we use 0.35 to determine whether we use Silhouette index as the criteria to determine the optimal number of clusters
    }
    
    # if no threshold for the height difference is provided
    if (missing(height.Ntimes) || is.null(height.Ntimes)) {
        height.Ntimes = 2  #by default, we select the first height which is (height.Ntimes) times larger than the immediate consecutive height
    }
    
    # just use simple criteria to determine whether they are feature vectors or
    # similarity matrix, and then we use different ways to measure the distance
    if (isSymmetric(mat)) {
        # symmmetric matrix
        d = as.dist(1 - mat)
        flag1 = 1
    } else {
        d = as.dist(1 - cor(t(mat)))
        flag1 = 0
    }
    
    h = hclust(d, method = hmethod)  #ward to ward.D
    
    
    # if N.cluster is given, we simply use the given N.cluster for hierarchical
    # clustering
#     if (!missing(N.cluster) && is.numeric(N.cluster)) {
    if (is.numeric(N.cluster)) {
        if (!is.numeric(N.cluster)) {
            stop("The given N.cluster is not a numeric!")
        }
        if (N.cluster%%1 != 0) {
            stop("The given N.cluster is not an integer!")
        }
        if (N.cluster < 2) {
            stop("The given N.cluster is less than 2, which is not suitable for clustering!")
        }
        
        v = cutree(h, k = N.cluster)  #for different numbers of clusters
        f = v  #the optimal clustering results
        sil = silhouette(v, d)
        msil = median(sil[, 3])
        ch0 = intCriteria(data.matrix(mat), as.integer(v), "Calinski_Harabasz")
        CHind = unlist(ch0, use.names = F)  #convert a list to a vector/value
        optN.cluster = N.cluster
    } else {
        # if missing, automatically determine the number of clusters
        
        my = mat
        nn = nrow(my)#number of data
        nc = minN.cluster:min(maxN.cluster, nrow(my) - 1)
#         cat("Testing numbers of clusters:", nc, "\n")
        # cat('trying cluster number as: ', nc)
        v = matrix(0, nrow = nrow(my), ncol = length(nc))  #for all different numbers of clusters
        msil = rep(0, length(nc))  #declare a vector of zeros
        # wss = rep(0, length(nc))#within-cluster sum of squares mdunn = rep(0, 39)#for
        # dunn index mdb = rep(0, 39)#for dunn index
        CHind = rep(0, length(nc))
        
        # print(paste('The height for the top 10 are: ', tail(h$height, n = 10), sep =
        # ''))
        
        my1 = as.matrix(my)  #convert to full matrix
        my = my1
        tt = numeric(length(nc))
        # cat(dim(my))
        for (i in 1:length(nc)) {
            # for(i in 1:1){ foreach(i=1:length(nc)) %dopar%{ cat('fast\n')
            
            v[, i] = cutree(h, k = nc[i])  #for different numbers of clusters
            
            sil = silhouette(v[, i], d)  #calculate the silhouette index
            
            # msil[i] = mean(sil[,3])#the mean value of the index
            msil[i] = median(sil[, 3])  #the mean value of the index
            # cat(Sys.time(), '\n') mdunn[i] = dunn(d, v[,i])
            
            # db = index.DB(d, cl = v[, i]) mdb[i] = db$DB
            
            # #within-cluster sum of squares spl <- split(d, v[,i]) wss[i] <- sum(sapply(spl,
            # wss))
            CHind[i] = get_CH(my, v[, i], disMethod = "1-corr")
            
#             if(nn>=100){
#                 tt[i] = length(which(table(v[, i]) < nn/100))#number of clusters whose number of cells is less than 50
#             }
            
            # ch0 = intCriteria(data.matrix(my),as.integer(v[,i]),'Calinski_Harabasz')
            # CHind[i] = unlist(ch0, use.names = F)#convert a list to a vector/value
            # cat(Sys.time(), '\n')
            
        }
        # cat('Should be not slow!')
        
        # print(msil)#the average sil index for all cases print(CHind)
        
        
        # print(wss) print(mdunn) print(mdb) oind = which.max(msil)#the index
        # corresponding to the max sil index
        tmp = which(msil == max(msil))  #in case there are more than one maximum
         
        if (length(tmp) > 1) {
            oind = tmp[ceiling(length(tmp)/2)]
        } else {
            oind = tmp
        }
        cat("The maximum Silhouette index is", max(msil), "\n")
        
        # if the maximum Silhouette index is smaller than the threshold, we use CH index
#         if (max(msil) > sil.thre){#pass through threshold
#             if(nrow(mat)>=100){
#             nt = which(tt == 0)#some clustering results having no cases of single-element as a cluster
#                 if(length(nt) != 0){
#                     x = rank(-msil[nt])
#                     oind = nt[x[1]]
#                 }else{
#                     cat("contain a cluster with very small number of data\n")
#                     y = rank(tt)
#                     oind = y[1]
#                 }
#             }
#            
# #             if(oind == 1){
# #                 x = rank(-msil)
# #                 a1 = max(table(v[, oind]))
# #                 a2 = min(table(v[, oind]))
# #                 if(nrow(mat)>=100 && a1/a2 > 10){#significant balance
# #                     oind = x[2]
# #                 }
# #             }
#         }else{
        if(max(msil) <= sil.thre){
            oind = which.max(CHind)
            if (oind == 1) {
                # if the maximum CH index with the minimum number of clusters, it's likely that
                # the CH index is not reliable either
                tmp = tail(h$height, n = 10)  #the height
                diftmp = diff(tmp)
                flag = diftmp > (height.Ntimes - 1) * tmp[1:(length(tmp) - 1)]  #require the height is more than (height.Ntimes) times of the immediate consecutive one
                
                if (any(flag)) {
                  # if any satifies the condition; make sure at least one satisfied
                  pind = which.max(flag)
                  opth = (tmp[pind] + tmp[pind + 1])/2  #the optimal height to cut
                  optv = cutree(h, h = opth)  #using the appropriate height to cut
                  oind = length(unique(optv)) - 1  #for consistency
                }
            }
            
            
            
            # difCH = diff(CHind) x1 = which(difCH<0) if(length(x1)>0){ oind = min(x1)#local
            # maximum }else{ oind = 1 }
            
        }
        # #if the silhouette index is not reliable, we use the gap statistics
        # if(max(msil)<=0.1){ gskmn <- clusGap(my, FUN = kmeans, nstart = 20, K.max = 15,
        # B = min(nrow(my), 50)) g = gskmn$Tab gap = g[, 'gap']#the gap info print(gap) #
        # oind = which.max(gap)#maximum gap sesim = g[, 'SE.sim']#standard error of the
        # gap print(sesim) oind = maxSE(gap, sesim)#maximize the gap with parsimony of
        # the model # if(oind >=floor(maxc*0.8)){#if the gap stastic keeps increasing
        # until very big number of cluster, we use the first SE-rule (Gap(k) >= Gap(k+1)
        # - SE) # sesim = g[, 'SE.sim']#standard error of the gap # print(sesim) # oind =
        # maxSE(gap, sesim)#maximize the gap with parsimony of the model # } }
        
        f = v[, oind]  #the optimal clustering results
        optN.cluster = length(unique(f))
        
    }
    
    
    hres = list()
    hres$f = f#optimal clustering results
    hres$v = v#different numbers of clustering results
    hres$maxsil = max(msil)
    hres$msil = msil
    hres$CHind = CHind
    hres$height = h$height
    hres$optN.cluster = optN.cluster
    
    return(hres)
}
