### Data loading stages

setStage("loadSample", "Load sample", "character", "cpSample")

setProtocol("xcms", "From file via xcms",
            representation(includeMSn = "logical"),
            function(object, includeMSn = FALSE) {
              cpSample(object, profstep = 0, includeMSn = includeMSn)
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

