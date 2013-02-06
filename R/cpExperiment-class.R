setClass("cpExperiment", representation(comps = "matrix"),
         contains = c("PipelineData", "xcmsSet"))
##' \code{cpExperiment} class exntends class \code{xcmsSet} and
##' \code{PipelineData}, used to construct \code{cpExperiment} object.
##'
##' It import LC/GC-MS files and passed to \code{xcmsSet} then attacha pipeline
##' to the returned object.
##' @title cosntructor for cpExperiment object
##' @aliases cpExperiment cpExperiment-class 
##' @param files path names of the NetCDF/mzXML files to read
##' @param snames sample names
##' @param sclass sample classes
##' @param phenoData sample names and classes
##' @param profmethod method to use for profile generation
##' @param profparam parameters to use for profile generation
##' @param nSlaves number of slaves/cores to be used for parallel peak
##' detection.  MPI is used if installed, otherwise the snow package is employed
##' for multicore support.
##' @param ... passed to function \code{xcmsSet} function.
##' @return a \code{cpExperiment} object.
##' @author Michael Lawrence, Tengfei Yin
##' @export cpExperiment cpExperiement-class
cpExperiment <- function(files = NULL, snames = NULL, sclass = NULL,
                         phenoData = NULL, profmethod = "bin",
                         profparam = list(), nSlaves = 0, ...)
{
  xset <- xcmsSet(files, sclass, snames, phenoData, profmethod, profparam,
                  nSlaves = nSlaves, ...)
  loadproto <- Protocol("loadSample", "xcms")
  pinfo <- profinfo(xset)
  args <- list("genProfile", "intbin")
  pinfo$integrate <- pinfo$method == "intlin"
  pinfo$method <- NULL
  pinfo$step <- 0 ## xcmsSet does not generate profile matrix -- up to findPeaks
  profproto <- do.call("Protocol", c(args, pinfo))
  peakargs <- list(...)
  method <- peakargs$method
  if (is.null(method))
    method <- "matchedFilter"
  peakargs$method <- NULL
  peakproto <- do.call("Protocol", c(list("findPeaks", method), peakargs))
  pipeline <- Pipeline(loadproto, profproto, peakproto)
  proto <- Protocol("loadExperiment", "xcms", pipeline = pipeline)
  new("cpExperiment", xset, pipeline = Pipeline(proto))  
}

#### Extract results as an ExpressionSet
setAs("cpExperiment", "ExpressionSet",
      function(from) {
        require(Biobase)
        result <- matrix(NA, nrow(from@groups), nrow(from@phenoData))
        comps <- from@comps
        lst <- by(comps,comps[,'group'],function(comps){
          rt <- median(comps[,'rt'])
        }
                  )
        rt <- unlist(lst)
        rt <- as.numeric(rt)/60
        comps <- from@comps[!is.na(from@comps[,"quantity"]),]
        comps <- from@comps[order(from@comps[,'group']),]
        groups <- as.numeric(factor(comps[,"group"]))
        result[(comps[,"sample"]-1)*nrow(result) + groups] <-
          log(comps[,"quantity"]+1)
        bad_groups <- apply(result, 1, function(x) all(is.na(x)))
        gp <- from@groups
        gp <- cbind(gp,rt)
        ## should I leave NA for missing values
        result <- result[!bad_groups,]
        result[is.na(result)] <- min(result[!is.na(result)])*0.9
        class <- from@phenoData$class
        features <- as(as.data.frame(gp[!bad_groups,]),
                       "AnnotatedDataFrame")
        new("ExpressionSet", exprs = result, featureData = features,
            phenoData = as(from@phenoData, "AnnotatedDataFrame"))
      })


setMethod("show", "cpExperiment", function(object) {
  callNextMethod()
  cat("\nVia a ")
  show(object@pipeline)
})
