### FIXME: make interactive

setMethod("image", "cpSample", function(x) {
  prof <- x@env$profile
  imaged <- data.frame(time = rep(x@scantime, each = nrow(prof)),
                       mz = rep(profMz(x), ncol(prof)),
                       intensity = as.vector(prof))
  if(min(prof)>=0) {minprof <- min(prof)}
  else(minprof <- 0)
  imaged <- imaged[prof>minprof,]
  plot_image(imaged,x)
})

setMethod("explore", c("xcmsRaw", "ProtoGenProfile"),
          function(object, protocol) 
          {
            p <- image(object)
            stack_plots("Profile Matrix Generating",list(p))
          })
