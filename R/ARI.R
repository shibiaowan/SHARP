#' Calculate the adjusted rand index
#'
#' This function is to calculate the performance of the algorithm by 5 metrics, including the Rand index, HA, MA, FM and Jaccard. 
#'
#' @param ground_true_clusters the ground-truth clusters
#' @param pred_clusters the predicted clusters
#' 
#' @examples
#' finalmetrics = ARI(gtc, finalrowColor)
#'
#' @export
ARI <- function(ground_true_clusters, pred_clusters){#w: the ground-truth clusters; rowColor: the predicted clusters

	truelabel = ground_true_clusters#get the label (categorical)
	
	rowColor = pred_clusters#the predicted clusters

	truec = truelabel#make a copy
	levels(truec) = c(1:length(levels(truec)))#convert the categorical to numeric

	truecl = as.numeric(as.character(truec))#convert to numeric vector
	
	#convert a (categorical) vector to a numeric vector 
	t1 = factor(rowColor)#rowColor is a vector like "red red purple brown ..."
	levels(t1) = c(1:length(levels(t1)))#convert the categorical to numeric
	t1 = as.numeric(as.character(t1))
	
	metrics = adjustedRand(truecl, t1)#the ARI performance metrics
	
	return(metrics)
}