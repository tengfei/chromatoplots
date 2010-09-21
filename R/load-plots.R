### TODO: marginal plots (TIC, total mass spectrum)
### TODO: querying image to extract chromatograms and spectral plots.
### TODO: image plot drawn as image (not as an R graphic, much too slow)
### TODO: support drawing two images at once (for comparison)
setMethod("explore", c("cpSample", protocolClass("loadSample")),
          function(object, protocol, ...)
          {
            p <- plotRaw(object, ...)
            stack_plots("Raw Data Input",list(p))
          })

setMethod("plotRaw", "cpSample", function(object, ...) {
  p <- plot_image(as.data.frame(rawMat(object, ...)))
  p
  
})
