setClass("cpSample", contains = c("PipelineData", "xcmsRaw"))

cpSample <- function(filename, profstep = 1, profmethod = "intlin",
                     profparam = list(), includeMSn = FALSE)
{
  xraw <- xcmsRaw(filename, profstep, profmethod, profparam, includeMSn)
  protocols <- list(Protocol("loadSample", "xcms", includeMSn = includeMSn))
  integrate <- profmethod == "intlin"
  args <- list("genProfile", "intbin", step = profstep, integrate = integrate)
  protocols <- c(protocols, do.call("Protocol", c(args, profparam)))
  new("cpSample", xraw, pipeline = new("Pipeline", protocols))
}

setMethod("show", "cpSample", function(object) {
  callNextMethod()
  cat("\nVia a ")
  show(object@pipeline)
})
