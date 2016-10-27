RocPlot <- function(data, legend = TRUE) {
  
  data <- data %>% dplyr::filter(!is.na(null))
  methods <- unique(data$method)
  n.methods <- length(methods) 
  title <- data$name[1]
  
  null <- data %>% dplyr::filter(method == methods[1]) %>% .$null
  alternative <- data %>% dplyr::filter(method == methods[1]) %>% .$alternative
  n.test <- length(null)
  pred <- prediction(c(null, alternative), c(rep(1, n.test), rep(0, n.test)))
  perf <- performance(pred, "tpr", "fpr")
  plot(perf, xlim = c(0, 1), ylim = c(0, 1), colorize = FALSE, lwd = 4, lty = 1, 
       main = title, cex.main = 2.2, cex.axis = 2, cex.lab = 1.8, 
       ylab = "TPR", xlab = "FPR")
  
  for (i in 2:n.methods) {
    null <- data %>% dplyr::filter(method == methods[i]) %>% .$null
    alternative <- data %>% 
      dplyr::filter(method == methods[i]) %>% .$alternative
    n.test <- length(null)
    pred <- prediction(c(null, alternative), c(rep(1, n.test), rep(0, n.test)))
    perf <- performance(pred, "tpr", "fpr")
    plot(perf, lwd = 4, lty = i, add = TRUE, col = i)  
  }
  
  if (legend) {
    legend("bottomright", legend = methods, 
           col = seq(1, n.methods), lwd = 4, lty = seq(1, n.methods), cex = 2)
  }
}


DistancePlot <- function(data) {
  return(ggplot(data, aes(method, distance)) + 
           geom_boxplot(aes(fill = method)) + facet_grid(.~ name))
}

Multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  library(grid)
  # From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/ 
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

