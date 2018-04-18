#' hierarchical clustering with number of clusters determined automatically
#'
#' This function is to do hierarchical clustering, with the number of clusters determined by a strategy combining Silhouette index, CH index and the heights after hierarchical clustering.
#'
#' @param E the feature matrix whose columns represent the features and whose rows represent data/cells.
#' @param tag the tag for each independent random projection
#' @param colorL the color set for different clusters
#'
#' @examples
#' rowColor= getrowColor(E1, tag, colorL)
#'
#' @import cluster
#'
#' @import clues
#'
#' @import clusterCrit
#' @export
getrowColor <- function(Emat){#hierarchical clustering
	
        my = Emat
	my = my[apply(my,1,function(x) sd(x)!=0),]
	my <- t(scale(t(my)))

        
        d=as.dist(1-cor(t(my)))
  
        hres = gethclust(d, my)
        f = hres$f
  
        nf = as.character(f)


        unf = unique(nf)
        N.cluster = length(unf)
	print(paste("The optimal number of clusters for individual RP is: ", N.cluster, sep = ""))
	
	
	rowColor=vector(mode="character",length=nrow(my))

	for(j in 1:N.cluster){
# 		sub.tf=cutree(h,k=N.cluster)==j
		sub.tf=nf==unf[j]
		if(j>length(colorL)){j = j %% length(colorL); if(j == 0){j = length(colorL)}}
		rowColor[sub.tf]=colorL[j]
# 		if(grepl("_", colorL[j])){#whether the color word has "_" or not
#                     rowColor2[sub.tf] = unlist(strsplit(colorL[j], "_"))[1]#this is to remove _1, _2 to use the color words
#                 }else{
#                     rowColor2[sub.tf] = rowColor[sub.tf]
#                 }
# 		clustername = paste(outdir, "Cluster_",str_replace(fname,".txt",""), tag, ".",colorL[j],".",j,".txt",sep="")
# 		print(clustername)
# 		write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
# 		#write.table(names(sub.tf[sub.tf==TRUE]),clustername,col.names=F,row.names=F)
# 		cluster[sub.tf]=j
		
	}
	
# 	f_seind = f[seind]#the corresponding clustering results
	
# 	write.table(rowColor, paste(outdir, "AllClusters_rowColor",str_replace(fname,".txt",""), tag, "_hclust.txt", sep = ""), quote = FALSE, sep = "\n", row.names = F, col.names = F)#save the clustering colors
	
	
# 	cluster <- as.matrix(cluster)
# 	rownames(cluster)<-rownames(my)
# 	imgDat = t(my[h$labels[h$order],])
# 
# 	if (flagArrangeColumn ==1) imgDat = imgDat[h2$labels[h2$order],]	# cluster column as well
# 	#P =      mlabel[h$labels[h$order],]
# 	imgDat[imgDat>maxval]=maxval
# 	imgDat[imgDat<minval]=minval
# 
# 	#color code
# 	#par(mar=c(10,33,5,0))
# 	mycol = colorL[1:N.cluster]
# 	cluster = cluster[h$labels[h$order],]
	#image(t(as.matrix(as.numeric(cluster))),col=mycol,  axes=FALSE)
	#box()
	#par(mar=c(10,1,5,1))
					
# 	#-------------------it has used the rowColor here in RowSideColors-----------------------------				
# 	heatmap.2(as.matrix(my), Rowv=dend, Colv=F, scale="row",trace="none",dendrogram="row",labRow=F,
# 		 col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=0.7, RowSideColors=rowColor2,
# 		 cexCol=1.5,margins=c(15,5) )#here RowSideColors=rowColor2 instead of rowColor

# 	dev.off()

# 	pngname2=paste(outdir,fname,".Dendro2.png",sep="")
# 	png(pngname2,width=1000,height=1200)
# 	    heatmap.2(as.matrix(my), Rowv=dend, Colv=F, scale="row",trace="none",dendrogram="row",labRow=F,
# 		 col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=0.7, #RowSideColors=F,
# 		 cexCol=1.5,margins=c(15,5) )
# 	dev.off()
	
# 	sil = silhouette (rowColor,d)
# 	plot(sil)
	
	res = list()
	res$rowColor = rowColor
	res$maxsil = hres$maxsil
# 	return(rowColor)
	return(res)
}

library(clues)
gethclust<- function(d, my){
    
	
        h=hclust(d, method="ward.D")#ward to ward.D
		
	#determine the optimal number of clusters	
	nc = 2:min(40, nrow(my)-1)
	v = matrix(0, nrow = nrow(my), ncol = length(nc))#for all different numbers of clusters
	msil = rep(0, length(nc))#declare a vector of zeros
# 	wss = rep(0, length(nc))#within-cluster sum of squares
# 	mdunn = rep(0, 39)#for dunn index
# 	mdb = rep(0, 39)#for dunn index
	CHind = rep(0, length(nc))
	
# 	print(paste("The height for the top 10 are: ", tail(h$height, n = 10), sep = ""))
	
	my1 = as.matrix(my)#convert to full matrix
	my = my1
	for(i in 1:length(nc)){
# 	foreach(i=1:length(nc)) %dopar%{
	  v[,i] = cutree(h, k = nc[i])#for different numbers of clusters
	  sil = silhouette(v[,i], d)#calculate the silhouette index 
# 	  msil[i] = mean(sil[,3])#the mean value of the index
	  msil[i] = median(sil[,3])#the mean value of the index
	  
# 	  mdunn[i] = dunn(d, v[,i])
	  
# 	  db = index.DB(d, cl = v[, i])
# 	  mdb[i] = db$DB

# 	  #within-cluster sum of squares
# 	  spl <- split(d, v[,i])
# 	  wss[i] <- sum(sapply(spl, wss))
	  CHind[i] = get_CH(my, v[,i], disMethod = "1-corr")
	}

# 	print(msil)#the average sil index for all cases
# 	print(CHind)
	
	
# 	print(wss)
# 	print(mdunn)
# 	print(mdb)
# 	oind = which.max(msil)#the index corresponding to the max sil index
	tmp = which(msil == max(msil))#in case there are more than one maximum
	if(length(tmp)>1){oind = tmp[ceiling(length(tmp)/2)]}else{oind = tmp}
	
	if(max(msil)<=0.35){
	  oind = which.max(CHind)
# 	  if(oind ==1){#it's likely that the CH index is not reliable either
# 	    tmp = tail(h$height, n =10)#the height
# 	    diftmp = diff(tmp)
# 	    flag = diftmp > tmp[1:(length(tmp)-1)]#require the height is more than 2 times of the immediate consecutive one
# 	    
# 	    if(any(flag)){#if any satifies the condition; make sure at least one satisfied
# 	      pind = which.max(flag)
# 	      opth = (tmp[pind] + tmp[pind+1])/2#the optimal height to cut
# 	      optv = cutree(h, h = opth)#using the appropriate height to cut
# 	      oind = length(unique(optv)) - 1#for consistency
# 	    }
# 	  }



# 	  difCH = diff(CHind)
# 	  x1 = which(difCH<0)
# 	  if(length(x1)>0){
# 	    oind = min(x1)#local maximum
# 	  }else{
# 	    oind = 1
# 	  }
	  
	}
# 	#if the silhouette index is not reliable, we use the gap statistics
# 	if(max(msil)<=0.1){
# 	  gskmn <- clusGap(my, FUN = kmeans, nstart = 20, K.max = 15, B = min(nrow(my), 50))
# 	  g = gskmn$Tab
# 	  gap = g[, "gap"]#the gap info
# 	  print(gap)
#     
# # 	  oind = which.max(gap)#maximum gap
# 	  
# 	  sesim = g[, "SE.sim"]#standard error of the gap
# 	  print(sesim)
# 	  oind = maxSE(gap, sesim)#maximize the gap with parsimony of the model
#     
# # 	  if(oind >=floor(maxc*0.8)){#if the gap stastic keeps increasing until very big number of cluster, we use the first SE-rule (Gap(k) >= Gap(k+1) - SE)
# # 	    sesim = g[, "SE.sim"]#standard error of the gap
# # 	    print(sesim)
# # 	    oind = maxSE(gap, sesim)#maximize the gap with parsimony of the model
# # 	  }
# 	}
	
# 	oind = 10
	f = v[,oind]#the optimal clustering results
    
        hres = list()
        hres$f = f
        hres$maxsil = max(msil)
        hres$msil = msil
        hres$CHind = CHind
        hres$height = h$height
        hres$oind = oind
        
        return(hres)
}