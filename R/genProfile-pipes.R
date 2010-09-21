## Link xcms profile generation into the pipeline

setStage("genProfile", "Generate profile matrix", "cpSample")

setProtocol("intbin", "Interpolated bins",
            representation(step = "numeric", baselevel = "numeric",
                           basespace = "numeric", integrate = "logical"),
            function(data, step = 1, baselevel = min(data@env$intensity),
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
