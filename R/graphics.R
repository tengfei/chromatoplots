### GENERAL/SHARED GRAPHICS UTILITIES
assignments_to_edges <- function(assignments) {
  do.call("rbind", lapply(unique(assignments[!is.na(assignments)]), 
                          function(assignment) 
                          {
                            g <- which(assignments == assignment)
                            cbind(g[-length(g)], g[-1])
                          }))
}

plot_image <- function(data,x)
{
  xlim <- range(x@scantime)
  ylim <- x@mzrange
  zlim <- c(0,log(max(x@env$intensity)))
  p <- ggplot(data, aes(x = time, y = mz, fill = log(intensity))) + geom_tile()
  p <- p + scale_fill_gradient(limits=zlim,low = "yellow", high = "black")+
    opts(legend.position="none")
  p <- p+scale_x_continuous(limits=xlim)+scale_y_continuous(limits=ylim)
  p
}

stack_plots <- function(title, plots, weights = rep(1, length(plots)))
{
  win <- gtkWindow(show = FALSE)
  win$setTitle(paste("Chromatoplots [Stage: ", title, "]", sep=""))
  win$setDefaultSize(800, 600)
  win$show()
  vbox <- gtkVBox()
  win$add(vbox)
  sizes <- weights / sum(weights) * win$getDefaultSize()$height
  for (i in seq_along(plots)) {
    plot <- plots[[i]]
    if (inherits(plot, "ggplot"))
      wid <- gtkDrawingArea()
    else # ggobi plot
      wid <- plot
    wid$showAll()
    vbox$add(wid)
    wid$setSizeRequest(-1, sizes[i])
    if (inherits(plot, "ggplot")) {
      asCairoDevice(wid)
      print(plot)
    }
  }
}
