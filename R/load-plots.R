### TODO: marginal plots (TIC, total mass spectrum)
### TODO: querying image to extract chromatograms and spectral plots.
### TODO: image plot drawn as image (not as an R graphic, much too slow)
### TODO: support drawing two images at once (for comparison)
setMethod("explore", c("cpSample", protocolClass("loadSample")),
          function(object, protocol,...)
          {
            p <- plotRaw(object, ...)
            stack_plots("Raw Data Input",list(p))
          })

setMethod("plotRaw", "cpSample", function(object,...) {
  xlim <- range(object@scantime)
  ylim <- range(object@env$mz)
  zlim <- c(0,log(max(object@env$intensity)))
  dfr <- as.data.frame(rawMat(object))
  int.cut <- quantile(dfr$intensity,0.95)
  dfr2 <- dfr[dfr$intensity>int.cut,]
  p <- ggplot(dfr2, aes(x = time, y = mz, fill = log(intensity))) + geom_tile()
  p <- p + scale_fill_gradient(limits=zlim,low = "yellow", high = "black")+
    opts(legend.position="none")
  p <- p+scale_x_continuous(limits=xlim)+scale_y_continuous(limits=ylim)
  p
})

