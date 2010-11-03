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


