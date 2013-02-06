setMethod("explore", c("cpSample", protocolClass("loadSample")),
          function(object, protocol,...)
          {
            p <- plotRaw(object, ...)
            stack_plots("Raw Data Input",list(p))
          })

plotRaw <- function(object, legend = TRUE) {
  ## xlim <- range(object@scantime)
  ## ylim <- range(object@env$mz)
  ## zlim <- c(0,log(max(object@env$intensity)))
  dfr <- as.data.frame(xcms:::rawMat(object))
  p <- ggplot(dfr, aes(x = time, y = mz, fill = log(intensity))) + geom_raster()
  ## int.cut <- quantile(dfr$intensity,0.95)
  ## dfr2 <- dfr[dfr$intensity>int.cut,]
  ## p <- ggplot(dfr2, aes(x = time, y = mz, fill = log(intensity))) + geom_raster()
  ## p <- p + scale_fill_gradient(limits=zlim,low = "yellow", high = "black")
  p <- p + scale_fill_gradient(low = "yellow", high = "black")
  if(!legend)
    p <- p + theme(legend.position="none")
  ## p <- p+scale_x_continuous(limits=xlim)+scale_y_continuous(limits=ylim)
  p <- p + theme_bw()
  message("Making rastered graphics ...")
  p
}

