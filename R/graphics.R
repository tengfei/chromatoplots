setGeneric('cplot',function(obj,protocol,...) standardGeneric('cplot'))

setMethod('cplot',c('PipelineData','missing'),
          function(obj,protocol,...){
            proto <- NULL
            if(length(obj@pipeline@.Data))
              proto <- tail(obj@pipeline,1)[[1]]
            cplot(obj,proto,...)
          })


### GENERAL/SHARED GRAPHICS UTILITIES
assignments_to_edges <- function(assignments) {
  do.call("rbind", lapply(unique(assignments[!is.na(assignments)]), 
                          function(assignment) 
                          {
                            g <- which(assignments == assignment)
                            cbind(g[-length(g)], g[-1])
                          }))
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
    vbox$add(wid)
    wid$showAll()
    wid$setSizeRequest(-1, sizes[i])
    if (inherits(plot, "ggplot")) {
      asCairoDevice(wid)
      print(plot)
    }
  }
}
