##' load samples
##'
##' load LC/GC-MS sample
##' @title load LC/GC-MS sample
##' @aliases loadSample-method
##' @author Michael Lawerence
setStage("loadSample", "Load sample", "character", "cpSample")

setProtocol("xcms", "From file via xcms",
            representation(includeMSn = "logical"),
            function(object, profstep = 1, profmethod = "intlin",
                     profparam = list(), includeMSn = FALSE,
                     mslevel = NULL, scanrange = NULL) {
              cpSample(object, profstep = profstep,
                       profmethod = profmethod,
                       profparam = profparam,
                       includeMSn = includeMSn,
                       mslevel = mslevel,
                       scanrange = scanrange)
            }, "loadSample")


setStage("loadExperiment", "Load experiment", "character", "cpExperiment")

setClassUnion("data.frameORNULL", c("data.frame", "NULL"))

setProtocol("xcms", "From file(s) via xcms",
            representation(pipeline = "Pipeline",
                           phenoData = "data.frameORNULL"),
            function(object, pipeline = LoadPeaksPipeline(), phenoData = NULL)
            { ### NOTE: pipeline does NOT control sample loading here
              peakpipe <- pipeline(pipeline, "cpSample", "cpPeaks")
              if (!length(peakpipe))
                stop("no pipe from cpSample to cpPeaks in 'pipeline'")
              cpExperiment(object, phenoData = phenoData,
                           method = "pipeline", pipeline = peakpipe)
            }, "loadExperiment")


setGeneric("processRawsProto",function(object,...)
           standardGeneric("processRawsProto"))

setMethod("processRawsProto","cpExperiment",function(object,
                                                     role="genProfile",
                                                     method=defaultMethod(role)){
  pc <- object@pipeline[[1]]@pipeline
  if(is.call(pc)) pps <- eval(pc)
  protos <- findProtocols(pps,role,method)
  if(!length(protos))
    NULL
  else{
    npps <- do.call("Pipeline",lapply(1:protos[1],function(i) pps[[i]]))
    return(npps)
  }
  
})

setMethod("processRawsProto","cpSample",function(object,
                                                 role="genProfile",
                                                 method=defaultMethod(role)){
  pc <- object@pipeline
  if(is.call(pc)) pps <- eval(pc)
  protos <- findProtocols(pps,role,method)
  if(!length(protos))
    NULL
  else{
    npps <- do.call("Pipeline",lapply(1:protos[1],function(i) pps[[i]]))
    return(npps)
  }
  
})


