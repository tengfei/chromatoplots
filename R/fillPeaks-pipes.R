setStage("fillPeaks", "Impute Missing Peaks", "cpExperiment")
setProtocol("extract", "Extract from raw data",
            fun = getMethodLocal("fillPeaks", "xcmsSet"),
            parent = "fillPeaks")
