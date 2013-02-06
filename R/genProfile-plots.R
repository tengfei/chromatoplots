### FIXME: make interactive
setMethod("image", "cpSample", function(x) {
  prof <- x@env$profile
  imaged <- data.frame(time = rep(x@scantime, each = nrow(prof)),
                       mz = rep(profMz(x), ncol(prof)),
                       intensity = as.vector(prof))
  if(min(prof)>=0) {
    minprof <- min(prof)
  }else{
    minprof <- 0
  }
  imaged <- imaged[prof>minprof,]
  ## imaged <- imaged[prof > 0,]

  plot_image(imaged,x)
})

setMethod("explore", c("xcmsRaw", "ProtoGenProfile"),
          function(object, protocol, ...) 
          {
            p <- image(object, ...)
            stack_plots("Profile Matrix Generating",list(p))
          })

plot_image <- function(data,x, legend = TRUE)
{
  xlim <- range(x@scantime)
  ylim <- x@mzrange
  zlim <- c(0,log(max(x@env$intensity)))
  p <- ggplot(data, aes(x = time, y = mz, fill = log(intensity))) + geom_raster()
  ## p <- p + scale_fill_gradient(limits=zlim,low = "yellow", high = "black")
  p <- p + scale_fill_gradient(low = "yellow", high = "black")
  if(!legend)
    p <-  p + theme(legend.position="none")
  p <- p + theme_bw()
  message("Making rastered graphics ...")  
  ## p <- p+scale_x_continuous(limits=xlim)+scale_y_continuous(limits=ylim)
  p
}


