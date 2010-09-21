setClass("cpPeaks", contains = c("PipelineData", "xcmsPeaks"))

setMethod("show", "cpPeaks", function(object) {
  callNextMethod()
  cat("\nVia a ")
  show(object@pipeline)
})
