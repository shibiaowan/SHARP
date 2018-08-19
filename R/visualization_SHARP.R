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
#' @param width, height the width and height of the outpuf figure.
#'
#' @examples
#' y = SHARP(scExp)
#' visualization_SHARP(y)
#'
#' @import Rtsne
#'
#' @import RColorBrewer
#'
#' @export

visualization_SHARP<- function(y, label, w, filename, filetype, width = 900, height = 900){

    # timing
    start_time <- Sys.time()  #we exclude the time for loading the input matrix
    cat("Start visualization...\n")
   

#     set.seed(10); x1 = jitter(y$x0, amount= 0)#to avoid overlapping
    
#     x1 = jitter(y$viE, amount= 0)#to avoid overlapping
    w1 = dim(y$x0)[2]
    w2 = dim(y$viE)[2]
    
    if(missing(w)){#the relative weight of x0 over ensemble RP-based matrices
        w = 2
    }
    x1 = cbind(w*scale(y$x0), scale(y$viE))
#     x1 = scale(y$viE)
#     x1 = cbind(y$x0, y$viE)
    #whether PCA is required
    if (dim(x1)[2] < 50){
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
    rtsne_out <- Rtsne(as.matrix(x1), check_duplicates = FALSE, pca = flag)
    cat("Project to 2-D space...\n")
    file_plot <- filename
    
    if (filetype == "pdf"){
        pdf(file_plot)
    }else if (filetype == "png"){
        png(file_plot, width = width, height = height)
    }
    
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
        palette(brewer.pal(n = uc, name = "Set1"))
        
        plot(rtsne_out$Y, asp = 1, pch = 20, col = label, cex = 0.75, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, xlab = "SHARP Dim-1", ylab = "SHARP Dim-2", main = tt)
    }else{
        plot(rtsne_out$Y, asp = 1, pch = 20, cex = 0.75, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5, xlab = "SHARP Dim-1", ylab = "SHARP Dim-2", main = tt)
    }
    
    dev.off()
    
    end_time <- Sys.time()
    
    t <- difftime(end_time, start_time, units = "mins")  #difference time in minutes
    cat("Running time for visualization:", t ,"minutes\n")
    cat("-----------------------------------------------------------------------\n")
    cat("Results saved into", file_plot, "\n")
}
