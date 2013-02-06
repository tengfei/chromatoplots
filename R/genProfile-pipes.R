## Link xcms profile generation into the pipeline
setStage("genProfile", "Generate profile matrix", "cpSample")
##' profile matrix method: intbin or binlinbase
##'
##' This stage generage profile matrix.
##' 
##' @title Profeile generation stage
##' @param data cpSample object.
##' @param step numeric value for m/z step.
##' @param baselevel the intensity level to fill in for empty parts of the
##' spectra. It defaluts to half of the minimum intensity.
##' @param basespace the m/z length after which the signal will drop to the base
##' level.
##' @param integrate if \code{FALSE} by default, use binlinbase profMethod, else
##' use intlin.
##' @param ... 
##' @return 
##' @author Michael Lawrence, Tengfei Yin
setProtocol("intbin", "Interpolated bins",
            representation(step = "numeric", baselevel = "numeric",
                           basespace = "numeric", integrate = "logical"),
            function(data, step = 1, baselevel = NULL,
                     basespace = 0, integrate = FALSE, ...)
            {
              if (integrate)
                method <- "intlin"
              else {
                method <- "binlinbase"
                data@profparam <- list(baselevel = baselevel,
                                       basespace = basespace)
              } 
              profMethod(data) <- method
              data@env <- copyEnv(data@env)
              profStep(data) <- step
              data
            }, "genProfile")
