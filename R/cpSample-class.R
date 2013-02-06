setClass("cpSample", contains = c("PipelineData", "xcmsRaw"))
## getting profile matrix too
##' \code{cpSample} is a constructor for creating \code{cpSample} object, it
##' imports single LC/GC-MS file.
##'
##' \code{cpSample} extends \code{xcmsRaw} class and \code{PipelineData}
##' class. So it's constrcutor looks identical to \code{xcmsRaw} constructor in
##' xcms package.
##' @title constructor for cpSample object
##' @aliases cpSample-class cpSample show,cpSample-method
##' @param filename path name of the NetCDF or mzXML file to read
##' @param profstep step size (in m/z) to use for profile generation
##' @param profmethod method to use for profile generation
##' @param profparam extra parameters to use for profile generation
##' @param includeMSn only for XML file formats: also read MS$^n$ (Tandem-MS of
##' Ion-/Orbi- Trap spectra)
##' @param mslevel move data from mslevel into normal MS1 slots, e.g. for peak
##' picking and visualisation
##' @param scanrange scan range to read
##' @return \code{cpSample} object.
##' @author Michael Lawrence, Tengfei Yin
##' @export cpSample cpSample-class
cpSample <- function(filename, profstep = 1, profmethod = "intlin",
                     profparam = list(), includeMSn = FALSE,
                     mslevel = NULL, scanrange = NULL)
{
  xraw <- xcmsRaw(filename, profstep, profmethod, profparam, includeMSn,
                  mslevel, scanrange)
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
