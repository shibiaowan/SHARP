#' SHARP visualization
#'
#' This function is to visualize the data in 2D scatter plots. For visualization, SHARP uses a weighted combination of the dimension-reduced feature matrix and the output location matrix produced by the clustering results of SHARP. The output scatter plot will be saved in a figure file (by default, the file is named as "vi_SHARP.pdf" or "vi_SHARP.png") 
#'
#' @param y the results after running the SHARP function
#'
#' @param label the reference (or pre-defined) clustering label for the data; if not given, the resulting plot will be drawn in black/white color.
#' 
#' @param w the relative weight of clustering-result based matrix over ensemble random-projection based matrix; by default, w = 2.
#'
#' @param filename the name of the file to which the output figure will be saved; if not given, the output figure will be saved in "vi_SHARP.pdf"
#'
#' @param filetype the type of the output file. Suggested file types are PDF or PNG, while other common types (e.g., JPEG or TIFF) are also acceptable. If not given, the type will be determined as follows: when the number of single cells is less than 5000, the type will be PDF; otherwise, it will be PNG.
#'
#' @param width/height the width/height of the outpuf figure.
#'
#' @examples
#' y = SHARP(scExp)
#' visualization_SHARP(y)
#'
#' @import Rtsne
#'
#' @import RColorBrewer
#'
#' @import parallel
#'
#' @import ggplot2
#'
#' @export

visualization_SHARP<- function(y, label, w, filename, filetype, n.cores, legendtitle = "Cell Type", width = 9.5, height = 8.5, res = 400, ...){

    # timing
    start_time <- Sys.time()  #we exclude the time for loading the input matrix
    if (missing(n.cores)) {
        # number of cores to be used, the default is to use all but one cores
        n.cores = detectCores() - 1
    }
    
    cat("Start visualization...\n")
   

#     set.seed(10); x1 = jitter(y$x0, amount= 0)#to avoid overlapping
    
#     x1 = jitter(y$viE, amount= 0)#to avoid overlapping
    w1 = dim(y$x0)[2]
    w2 = dim(y$viE)[2]
    
    if(missing(w)){#the relative weight of x0 over ensemble RP-based matrices
        w = 2
    }
    
    if(w >= 100){#if the weight is very large, only the clustering results related matrix is used
        x1 = as.matrix(y$x0)
        x1 = jitter(x1, amount = 0)
    }else if(w <= 0.01){#if the weight is very small, only the RP-projected matrix is used
        x1 = as.matrix(y$viE)
    }else{
        x1 = as.matrix(cbind(w*scale(y$x0), scale(y$viE)))
    }
    
#     x1 = scale(y$viE)
#     x1 = cbind(y$x0, y$viE)
    #whether PCA is required
    if (dim(x1)[2] <= 50){
        flag = FALSE
    }else{
        flag = TRUE
    }
    
    #if the file type is not given
    if(missing(filetype)){
        if(dim(x1)[1] < 5000){
            filetype = "pdf"
        }else{
            filetype = "png"
        }
    }
   
    #if the file name is not given
    if(missing(filename)){
        filename = paste("vi_SHARP.", filetype, sep = "")
    }
    
    set.seed(10) # Sets seed for reproducibility
    
#     if(dim(x1)[1] > 1e5){
#         rM = ranM(t(x1), 30, seedn = 61)
#         x1 = 1/sqrt(30) * t(rM) %*% t(x1)
#         x1 = t(x1)
#     }
    cat("Project to 2-D space...\n")
#     rtsne_out <- Rtsne(x1, check_duplicates = FALSE, pca = flag, num_threads = n.cores, verbose = TRUE, ...)
    rtsne_out <- Rtsne(x1, check_duplicates = FALSE, pca = flag, num_threads = n.cores, ...)
        file_plot <- filename
    
   
#     par(xpd = T, mar = par()$mar + c(0,0,0,7))
   
    tt = "2D SHARP Visualization"
    
#     allcol =  c("purple",  "pink",  "black",  "orange", "turquoise", "yellow", "beige", "gray", 
#          "coral",    "khaki",  "violet", "magenta",
#          "salmon", "goldenrod", "orchid", "seagreen", "slategray", "darkred", 
#         "darkblue", "darkcyan", "darkgreen", "darkgray", "darkkhaki", "darkorange", 
#         "darkmagenta", "darkviolet", "darkturquoise", "darksalmon", "darkgoldenrod", 
#         "darkorchid", "darkseagreen", "darkslategray", "deeppink", "lightcoral", 
#         "lightcyan")

    cat("Draw the scatter plots...\n")
    if(!missing(label)){#if the reference clustering label is given
        uc = length(unique(label))
#         palette(brewer.pal(n = uc, name = "Set1"))
        
#         d0 = data.frame(as.character(label))
#         colnames(d0) = c("label")
#         xa = max(nchar(as.character(unique(label))))#the maximum number of characters
#         cat("Max nchar:", xa, "\n")
#         #only in this case, we will draw the legend
#         par(mar=c(5.1, 4.1, 4.1, 2.1 + max(6, xa/2)), xpd=TRUE)
#         plot(rtsne_out$Y, asp = 1, pch = 20, col = d0$label, cex = 0.75, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, xlab = "SHARP Dim-1", ylab = "SHARP Dim-2", main = tt)
# #         legend(0.03, 0.015, legend = levels(d0$label), col = 1:length(d0$label), pch = 20)
#         legend("topright", inset=c(min(-0.1*xa/3.5, -0.2),0), legend = levels(d0$label), col = 1:length(d0$label), pch = 20)
#         
        
        d0 = data.frame(rtsne_out$Y, as.character(label))
        colnames(d0) = c("x1", "x2", "label")
        
        allcol = c("black", "red",  "green", "blue", "cyan", "magenta", "yellow", "grey", "brown", "purple",  "orange", "turquoise",  "beige",  
         "coral",    "khaki",  "violet", "pink", "salmon", "goldenrod", "orchid", "seagreen", "slategray", "darkred",  "darkblue", "darkcyan", "darkgreen", "darkgray", "darkkhaki", "darkorange", "darkmagenta", "darkviolet", "darkturquoise", "darksalmon", "darkgoldenrod", "darkorchid", "darkseagreen", "darkslategray", "deeppink", "lightcoral", "lightcyan")
        nl = length(allcol)
        n0 = 1:uc %% length(allcol)
        n0[n0 == 0] = nl#replace those 0's'
        pcol = allcol[n0]#recursively use the colors
        colScale <- scale_colour_manual(name = "label",values = pcol)#assign colors by ourselves
        
        vplot = ggplot(d0, aes(x = x1, y = x2, colour = label, group = label)) +
        theme_bw(base_size = 14)+
        theme_classic() +
        geom_point(size = 1) + 
        theme(axis.title=element_text(face="bold",size="14"), axis.text.x = element_text(size="14", hjust = 0.5, face="bold", colour= "black"), axis.text.y = element_text(size="14", face="bold", colour= "black"), legend.text=element_text(size="14", face="bold"), legend.title=element_blank(), plot.title = element_text(hjust = 0.5, size="14", face="bold"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA)) + xlab("SHARP Dim-1") + ylab("SHARP Dim-2") + labs(group = legendtitle) + ggtitle(tt) + colScale
        
        #legend.title=element_text(size="14")#with legend title
        
#         vplot = vplot + scale_fill_discrete(name = legendtitle)
        #        labs(fill = legendtitle) #### +
        
#         print(vplot)
        if (filetype == "pdf"){
#             pdf(file_plot, width = width)
            ggsave(file_plot, vplot, device = filetype, width = width)
        }else if (filetype == "png"){
#             png(file_plot, width = width, height = height, units = "in", res = res, pointsize = 4)
            ggsave(file_plot, vplot, device = filetype, units = "in", dpi = res)
        }
        
    }else{
        if (filetype == "pdf"){
            pdf(file_plot, width = width)
        }else if (filetype == "png"){
            png(file_plot, width = width, height = height, units = "in", res = res)
        }
        
        plot(rtsne_out$Y, asp = 1, pch = 20, cex = 0.75, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, xlab = "SHARP Dim-1", ylab = "SHARP Dim-2", main = tt)
        
        dev.off()
    }
    
#     par(mar=c(5, 4, 4, 2) + 0.1)#the default axis
   
    
    end_time <- Sys.time()
    
    t <- difftime(end_time, start_time, units = "mins")  #difference time in minutes
    cat("Running time for visualization:", t ,"minutes\n")
    cat("-----------------------------------------------------------------------\n")
    cat("Results saved into", file_plot, "\n")
}
