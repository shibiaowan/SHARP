Run_Mtimes_SHARP <- function(scExp, ensize.K, reduced.ndim, partition.ncells, base.ncells, n.cores, Mtimes){
rm(list=ls())#remove all 

##########Input matrix files###########

# Files = c("Wang/Pancreas_Julia.published.dataTpms_selected.txt")
# Files = c("Enge/tpms4MATLAB_selected.txt")#you need to make transpose of it because it is feature * samples
# Files = c("Goolam/Goolam_expression_log10.txt")
# Files = c("Park/Park_CPM.rds")
#########################################
A = Mtimes
K = ensize.K

# eps = c(0.3, 0.25, 0.2, 0.15, 0.1, 0.05)
allresults <- vector("list", length =  length(K))#declare a matrix of lists

for(k in K){
#   allresults = list()
#   for(i in 1:length(eps)){
    info = vector("list", length = A)
    kname = paste("enSize_", k, sep = "")
#     dim = ceiling(log2(2282)/(eps[i]^2))
#     kname = paste("dim_", dim, sep = "")
    for(j in 1:A){#times of using SHARP
      print("=========================================================================", quote = FALSE)
      print(paste("SHARP: Ensemble Size-", k, ", Run-", j, sep = ""), quote = FALSE)
# 	enresults = getallresults_large(Files, K = k)#using the default Dim
	enresults = SHARP(scExp, k, reduced.ndim, partition.ncells, base.ncells, n.cores)
      jname = paste("Run_", j, sep = "")
      info[[jname]] = enresults
    }
    allresults[[kname]] = info
#   }
# #   dset = enresults$dataname#just use the last one to retrieve the dataset name
#   
}


    print("Done!", quote = FALSE)
    return(allresults)
}